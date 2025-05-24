use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;

use aho_corasick::AhoCorasick;
use clap::Parser;
use rayon::prelude::*;

/// CLI arguments
#[derive(Parser)]
#[command(author, version, about="ddRAD simulator in Rust (no plot)")]
struct Cli {
    /// TSV of enzymes: name <tab> motif^cut
    enzymes_tsv: PathBuf,
    /// TSV of combos: combo_name <tab> enzA,enzB
    combos_tsv:   PathBuf,
    /// Reference FASTA
    fasta:        PathBuf,
    /// Minimum fragment length (inclusive)
    #[arg(long)]
    min: usize,
    /// Maximum fragment length (inclusive)
    #[arg(long)]
    max: usize,
    /// Number of threads to use
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,
    /// Output TSV for combo totals
    #[arg(long, default_value = "ddrad_totals.tsv")]
    totals_out:  PathBuf,
    /// Output TSV for per-chromosome summary
    #[arg(long, default_value = "ddrad_summary.tsv")]
    summary_out: PathBuf,
}

fn iupac_map() -> HashMap<char, Vec<char>> {
    let mut m = HashMap::new();
    m.insert('A', vec!['A']);
    m.insert('C', vec!['C']);
    m.insert('G', vec!['G']);
    m.insert('T', vec!['T']);
    m.insert('R', vec!['A','G']);
    m.insert('Y', vec!['C','T']);
    m.insert('S', vec!['G','C']);
    m.insert('W', vec!['A','T']);
    m.insert('K', vec!['G','T']);
    m.insert('M', vec!['A','C']);
    m.insert('B', vec!['C','G','T']);
    m.insert('D', vec!['A','G','T']);
    m.insert('H', vec!['A','C','T']);
    m.insert('V', vec!['A','C','G']);
    m.insert('N', vec!['A','C','G','T']);
    m
}

fn expand_iupac(motif: &str, map: &HashMap<char, Vec<char>>) -> Vec<String> {
    // use the 'N' entry as our default pool
    let default_pool = map.get(&'N').unwrap();
    let pools: Vec<&Vec<char>> = motif
        .chars()
        .map(|c| map.get(&c).unwrap_or(default_pool))
        .collect();

    let mut res = vec![String::new()];
    for pool in pools {
        let mut next = Vec::new();
        for prefix in &res {
            for &b in pool {
                let mut s = prefix.clone();
                s.push(b);
                next.push(s);
            }
        }
        res = next;
    }
    res
}

fn revcomp(s: &str) -> String {
    s.chars()
     .rev()
     .map(|c| match c {
         'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A',
         _   => 'N',
     })
     .collect()
}

fn parse_enzymes(path: &PathBuf) -> io::Result<Vec<(String,String,usize)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut enzymes = Vec::new();
    for line in reader.lines() {
        let ln = line?;
        let ln = ln.trim();
        if ln.is_empty() || ln.starts_with('#') { continue }
        let parts: Vec<&str> = ln.split_whitespace().collect();
        if parts.len() < 2 || !parts[1].contains('^') {
            eprintln!("skipping malformed enzyme line: {}", ln);
            continue;
        }
        let name = parts[0].to_string();
        let raw  = parts[1];
        let cut_idx = raw.find('^').unwrap();
        let motif   = raw.replace('^', "");
        enzymes.push((name, motif, cut_idx));
    }
    Ok(enzymes)
}

fn parse_combos(
    path: &PathBuf,
    enz_names: &HashMap<String,usize>
) -> io::Result<Vec<(String,usize,usize)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut combos = Vec::new();
    for line in reader.lines() {
        let ln = line?;
        let ln = ln.trim();
        if ln.is_empty() || ln.starts_with('#') { continue }
        let parts: Vec<&str> = ln.split_whitespace().collect();
        if parts.len() < 2 || !parts[1].contains(',') {
            eprintln!("skipping malformed combo line: {}", ln);
            continue;
        }
        let combo = parts[0].to_string();
        let (a,b)  = {
            let mut sp = parts[1].splitn(2, ',');
            (sp.next().unwrap(), sp.next().unwrap())
        };
        if let (Some(&i1), Some(&i2)) = (enz_names.get(a), enz_names.get(b)) {
            combos.push((combo, i1, i2));
        } else {
            eprintln!("unknown enzyme in combo: {}", ln);
        }
    }
    Ok(combos)
}

fn read_fasta(path: &PathBuf) -> io::Result<Vec<(String,String)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut recs = Vec::new();
    let mut name = String::new();
    let mut seq  = String::new();
    for line in reader.lines() {
        let ln = line?;
        if ln.starts_with('>') {
            if !name.is_empty() {
                recs.push((name.clone(), seq.clone()));
            }
            name = ln[1..].split_whitespace().next().unwrap().to_string();
            seq.clear();
        } else {
            seq.push_str(ln.trim());
        }
    }
    if !name.is_empty() {
        recs.push((name, seq));
    }
    Ok(recs)
}

fn count_ddrad_fragments_between(
    cuts1: &Vec<usize>,
    cuts2: &Vec<usize>,
    min_len: usize,
    max_len: usize
) -> usize {
    if cuts1.is_empty() || cuts2.is_empty() {
        return 0;
    }
    let (a, b) = if cuts1[0] < cuts2[0] {
        (cuts1, cuts2)
    } else {
        (cuts2, cuts1)
    };
    let mut i = 0;
    let mut j = 0;
    let mut count = 0;
    while i < a.len() && j < b.len() {
        while i < a.len() && a[i] <= b[j] {
            i += 1;
        }
        if i >= a.len() {
            break;
        }
        let len1 = b[j] - a[i - 1];
        if len1 >= min_len && len1 <= max_len {
            count += 1;
        }

        while j < b.len() && b[j] <= a[i] {
            j += 1;
        }
        if j >= b.len() {
            break;
        }
        let len2 = a[i] - b[j - 1];
        if len2 >= min_len && len2 <= max_len {
            count += 1;
        }
    }
    count
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    // Load enzymes and index map
    let enzymes = parse_enzymes(&cli.enzymes_tsv)?;
    let mut enz_name_to_idx = HashMap::new();
    for (i, (n, _, _)) in enzymes.iter().enumerate() {
        enz_name_to_idx.insert(n.clone(), i);
    }

    // Load combos
    let combos = parse_combos(&cli.combos_tsv, &enz_name_to_idx)?;

    // Build Aho–Corasick automaton
    let iupac = iupac_map();
    let mut patterns = Vec::new();
    let mut infos    = Vec::new();
    for (idx, (_, motif, cut_idx)) in enzymes.iter().enumerate() {
        for pat in expand_iupac(motif, &iupac) {
            let pat_len = pat.len();
            patterns.push(pat.clone());
            infos.push((idx, *cut_idx, pat_len));
            let rc = revcomp(&pat);
            if rc != pat {
                patterns.push(rc.clone());
                infos.push((idx, pat_len - cut_idx, pat_len));
            }
        }
    }
    let ac = AhoCorasick::new(&patterns);

    // Read FASTA
    let records = read_fasta(&cli.fasta)?;

    // Parallel digest & collect (combo, chromosome, count)
    let summary: Vec<(String,String,usize)> = records
        .par_iter()
        .flat_map_iter(|(chrom, seq)| {
            let mut cuts: Vec<Vec<usize>> = vec![Vec::new(); enzymes.len()];
            for mat in ac.find_iter(seq) {
                let (e, cut_i, pat_len) = infos[mat.pattern()];
                let start = mat.end() - pat_len;
                cuts[e].push(start + cut_i);
            }
            for v in &mut cuts {
                v.sort_unstable();
            }
            combos.iter().map(move |(combo, i1, i2)| {
                let cnt = count_ddrad_fragments_between(&cuts[*i1], &cuts[*i2], cli.min, cli.max);
                (combo.clone(), chrom.clone(), cnt)
            })
        })
        .collect();

    // Compute totals
    let mut totals_map: HashMap<String, usize> = HashMap::new();
    for (combo, _, cnt) in &summary {
        *totals_map.entry(combo.clone()).or_default() += *cnt;
    }
    let mut combo_order: Vec<String> = totals_map.keys().cloned().collect();
    combo_order.sort_by_key(|c| usize::MAX - totals_map[c]);

    // Write combo totals
    let mut fo = File::create(&cli.totals_out)?;
    writeln!(fo, "combo\ttotal_count")?;
    for combo in &combo_order {
        writeln!(fo, "{}\t{}", combo, totals_map[combo])?;
    }
    eprintln!("Wrote combo totals to {:?}", cli.totals_out);

    // Group summary per chromosome
    let mut summary_map: HashMap<String, Vec<(String,usize)>> = HashMap::new();
    for (combo, chrom, cnt) in summary {
        summary_map.entry(combo).or_default().push((chrom, cnt));
    }

    // Write per-chromosome summary
    let mut fo2 = File::create(&cli.summary_out)?;
    writeln!(fo2, "combo\tchromosome\tcount")?;
    for combo in &combo_order {
        if let Some(v) = summary_map.get(combo) {
            let mut v = v.clone();
            v.sort_by(|a,b| a.0.cmp(&b.0));
            for (chrom, cnt) in v {
                writeln!(fo2, "{}\t{}\t{}", combo, chrom, cnt)?;
            }
        }
    }
    eprintln!("Wrote per-chromosome summary to {:?}", cli.summary_out);

    Ok(())
}
