use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;

use aho_corasick::AhoCorasick;
use clap::Parser;
use rayon::prelude::*;

#[derive(Parser)]
#[command(author, version, about="ddRAD simulator in Rust")]
struct Cli {
    enzymes_tsv: PathBuf,
    combos_tsv:  PathBuf,
    fasta:       PathBuf,
    #[arg(long)] min: usize,
    #[arg(long)] max: usize,
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,
    #[arg(long, default_value = "ddrad_totals.tsv")]
    totals_out: PathBuf,
    #[arg(long, default_value = "ddrad_summary.tsv")]
    summary_out: PathBuf,
}

fn iupac_map() -> HashMap<char, Vec<char>> {
    let mut m = HashMap::new();
    m.insert('A', vec!['A']); m.insert('C', vec!['C']);
    m.insert('G', vec!['G']); m.insert('T', vec!['T']);
    m.insert('R', vec!['A','G']); m.insert('Y', vec!['C','T']);
    m.insert('S', vec!['G','C']); m.insert('W', vec!['A','T']);
    m.insert('K', vec!['G','T']); m.insert('M', vec!['A','C']);
    m.insert('B', vec!['C','G','T']); m.insert('D', vec!['A','G','T']);
    m.insert('H', vec!['A','C','T']); m.insert('V', vec!['A','C','G']);
    m.insert('N', vec!['A','C','G','T']);
    m
}

fn expand_iupac(motif: &str, map: &HashMap<char, Vec<char>>) -> Vec<String> {
    let default_pool = &map[&'N'];
    let pools: Vec<&Vec<char>> = motif
        .chars()
        .map(|c| map.get(&c).unwrap_or(default_pool))
        .collect();

    let mut acc = vec![String::new()];
    for pool in pools {
        let mut next = Vec::new();
        for prefix in &acc {
            for &base in pool {
                let mut s = prefix.clone();
                s.push(base);
                next.push(s);
            }
        }
        acc = next;
    }
    acc
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
    let mut v = Vec::new();
    for line in reader.lines() {
        let ln = line?.trim().to_string();
        if ln.is_empty() || ln.starts_with('#') { continue; }
        let parts: Vec<&str> = ln.split_whitespace().collect();
        if parts.len()<2 || !parts[1].contains('^') {
            eprintln!("Skipping malformed enzyme: {}", ln);
            continue;
        }
        let cut_idx = parts[1].find('^').unwrap();
        let motif   = parts[1].replace('^', "");
        v.push((parts[0].to_string(), motif, cut_idx));
    }
    Ok(v)
}

fn parse_combos(
    path: &PathBuf,
    idx_map: &HashMap<String,usize>
) -> io::Result<Vec<(String,usize,usize)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut v = Vec::new();
    for line in reader.lines() {
        let ln = line?.trim().to_string();
        if ln.is_empty() || ln.starts_with('#') { continue; }
        let parts: Vec<&str> = ln.split_whitespace().collect();
        if parts.len()<2 || !parts[1].contains(',') {
            eprintln!("Skipping malformed combo: {}", ln);
            continue;
        }
        let mut sp = parts[1].splitn(2, ',');
        let a = sp.next().unwrap();
        let b = sp.next().unwrap();
        if let (Some(&i1), Some(&i2)) = (idx_map.get(a), idx_map.get(b)) {
            v.push((parts[0].to_string(), i1, i2));
        } else {
            eprintln!("Unknown enzyme in combo: {}", ln);
        }
    }
    Ok(v)
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
    if cuts1.is_empty() || cuts2.is_empty() { return 0; }
    let (a,b) = if cuts1[0]<cuts2[0] {(cuts1,cuts2)} else {(cuts2,cuts1)};
    let mut i = 0;
    let mut j = 0;
    let mut cnt = 0;
    while i < a.len() && j < b.len() {
        while i < a.len() && a[i] <= b[j] { i += 1; }
        if i >= a.len() { break; }
        let d1 = b[j] - a[i-1];
        if (min_len..=max_len).contains(&d1) { cnt += 1; }

        while j < b.len() && b[j] <= a[i] { j += 1; }
        if j >= b.len() { break; }
        let d2 = a[i] - b[j-1];
        if (min_len..=max_len).contains(&d2) { cnt += 1; }
    }
    cnt
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    let enzymes = parse_enzymes(&cli.enzymes_tsv)?;
    let mut idx_map = HashMap::new();
    for (i,(name,_,_)) in enzymes.iter().enumerate() {
        idx_map.insert(name.clone(), i);
    }

    let combos = parse_combos(&cli.combos_tsv, &idx_map)?;

    let iupac = iupac_map();
    let mut patterns = Vec::new();
    let mut infos    = Vec::new();
    for (ei,(_, motif, cut_idx)) in enzymes.iter().enumerate() {
        for pat in expand_iupac(motif, &iupac) {
            let pat_len = pat.len();
            patterns.push(pat.clone());
            infos.push((ei, *cut_idx, pat_len));
            let rc = revcomp(&pat);
            if rc != pat {
                patterns.push(rc.clone());
                infos.push((ei, pat_len - cut_idx, pat_len));
            }
        }
    }
    let ac = AhoCorasick::new(&patterns);

    let records = read_fasta(&cli.fasta)?;

    let summary: Vec<(String,String,usize)> = records
        .par_iter()
        .flat_map_iter(|(chrom, seq)| {
            let mut cuts: Vec<Vec<usize>> = vec![Vec::new(); enzymes.len()];
            for mat in ac.find_overlapping_iter(seq) {
                let (e, cut_idx, _) = infos[mat.pattern()];
                let pos = mat.start() + cut_idx;
                cuts[e].push(pos);
            }
            for v in &mut cuts { v.sort_unstable(); }
            combos.iter().map(move |(combo, i1, i2)| {
                let cnt = count_ddrad_fragments_between(
                    &cuts[*i1], &cuts[*i2], cli.min, cli.max
                );
                (combo.clone(), chrom.clone(), cnt)
            })
        })
        .collect();

    let mut totals: HashMap<String, usize> = HashMap::new();
    for (combo, _, cnt) in &summary {
        *totals.entry(combo.clone()).or_default() += *cnt;
    }
    let mut combo_order: Vec<String> = totals.keys().cloned().collect();
    combo_order.sort_by_key(|c| usize::MAX - totals[c]);

    let mut fo = File::create(&cli.totals_out)?;
    writeln!(fo, "combo\ttotal_count")?;
    for combo in &combo_order {
        writeln!(fo, "{}\t{}", combo, totals[combo])?;
    }
    eprintln!("Wrote combo totals ➜ {:?}", cli.totals_out);

    let mut summary_map: HashMap<String, Vec<(String,usize)>> = HashMap::new();
    for (combo, chrom, cnt) in summary {
        summary_map.entry(combo).or_default().push((chrom, cnt));
    }

    let mut fo2 = File::create(&cli.summary_out)?;
    writeln!(fo2, "combo\tchromosome\tcount")?;
    for combo in &combo_order {
        if let Some(vec) = summary_map.get(combo) {
            let mut vec = vec.clone();
            vec.sort_by(|a,b| a.0.cmp(&b.0));
            for (chrom, cnt) in vec {
                writeln!(fo2, "{}\t{}\t{}", combo, chrom, cnt)?;
            }
        }
    }
    eprintln!("Wrote per-chromosome summary ➜ {:?}", cli.summary_out);

    Ok(())
}
