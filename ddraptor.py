#!/usr/bin/env python3
import sys
import re
import argparse
import multiprocessing as mp
from functools import partial
from itertools import product

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData
import ahocorasick
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

def expand_iupac(motif):
    table = IUPACData.ambiguous_dna_values
    pools = [table.get(b, b) for b in motif]
    return [''.join(p) for p in product(*pools)]

def parse_enzymes(path):
    enzymes = []
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        parts = re.split(r'\s+', ln)
        if len(parts) < 2 or '^' not in parts[1]:
            print(f"Warning: skipping malformed enzyme line: {ln}", file=sys.stderr)
            continue
        name, raw = parts[0], parts[1]
        cut_idx = raw.index('^')
        motif   = raw.replace('^','')
        enzymes.append((name, motif, cut_idx))
    return enzymes

def parse_combos(path):
    combos = {}
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        parts = re.split(r'\s+', ln)
        if len(parts) < 2 or ',' not in parts[1]:
            print(f"Warning: skipping malformed combo line: {ln}", file=sys.stderr)
            continue
        combo, pair = parts[0], parts[1]
        a, b = pair.split(',',1)
        combos[combo] = (a, b)
    return combos

def build_automaton(enzymes):
    A = ahocorasick.Automaton()
    for idx, (name, motif_iupac, cut_idx) in enumerate(enzymes):
        for seq_pat in expand_iupac(motif_iupac):
            L = len(seq_pat)
            # forward
            A.add_word(seq_pat, (idx, cut_idx, L))
            # reverse‐complement
            rc = str(Seq(seq_pat).reverse_complement())
            if rc != seq_pat:
                rc_cut = L - cut_idx
                A.add_word(rc, (idx, rc_cut, L))
    A.make_automaton()
    return A

def ddpositions(list1, list2):
    if not list1 or not list2:
        return []
    if list1[0] < list2[0]:
        start, end = list1, list2
    else:
        start, end = list2, list1
    i = j = 0
    dd = []
    while i < len(start) and j < len(end):
        while i < len(start) and start[i] <= end[j]:
            i += 1
        if i == len(start):
            break
        dd.append(start[i-1])
        dd.append(end[j])
        while j < len(end) and end[j] <= start[i]:
            j += 1
        if j == len(end):
            break
        dd.append(end[j-1])
        dd.append(start[i])
    if len(dd) % 2 != 0:
        dd = dd[:-1]
    return dd

def count_ddrad_fragments(record, enzymes, automaton, combos, min_len, max_len):
    seq = str(record.seq)
    cuts = { name: [] for name,_,_ in enzymes }
    for end_pos, (idx, cut_idx, L) in automaton.iter(seq):
        start_pos = end_pos - L + 1
        cuts[enzymes[idx][0]].append(start_pos + cut_idx)
    for v in cuts.values():
        v.sort()

    results = []
    for combo, (e1,e2) in combos.items():
        d = ddpositions(cuts[e1], cuts[e2])
        cnt = 0
        for k in range(0, len(d), 2):
            length = d[k+1] - d[k]
            if min_len <= length <= max_len:
                cnt += 1
        results.append((combo, record.id, cnt))
    return results

def main():
    p = argparse.ArgumentParser(
        description="ddRAD simulator with customizable outputs"
    )
    p.add_argument('enzymes_tsv', help="TSV: enzyme_name <tab> motif^cut")
    p.add_argument('combos_tsv',   help="TSV: combo_name <tab> enzA,enzB")
    p.add_argument('fasta',        help="Reference FASTA")
    p.add_argument('--min',  type=int, required=True,  help="Min fragment length")
    p.add_argument('--max',  type=int, required=True,  help="Max fragment length")
    p.add_argument('--processes', type=int, default=mp.cpu_count(),
                   help="Worker processes")
    p.add_argument('--totals-out',   default='ddrad_totals.tsv',
                   help="Combo totals TSV path")
    p.add_argument('--summary-out',  default='ddrad_summary.tsv',
                   help="Per-chromosome summary TSV path")
    p.add_argument('--heatmap-out',  default='ddrad_heatmap.png',
                   help="Heatmap PNG path")
    args = p.parse_args()

    enzymes = parse_enzymes(args.enzymes_tsv)
    combos  = parse_combos(args.combos_tsv)
    autom   = build_automaton(enzymes)

    worker = partial(count_ddrad_fragments,
                     enzymes=enzymes,
                     automaton=autom,
                     combos=combos,
                     min_len=args.min,
                     max_len=args.max)

    records = list(SeqIO.parse(args.fasta, "fasta"))
    pool    = mp.Pool(args.processes)
    all_rows = []
    for rec_res in tqdm(pool.imap(worker, records),
                        total=len(records), desc="Digesting"):
        all_rows.extend(rec_res)
    pool.close(); pool.join()

    df = pd.DataFrame(all_rows, columns=['combo','chromosome','count'])

    # 1) Combo totals
    totals = df.groupby('combo')['count'].sum().sort_values(ascending=False)
    totals_df = totals.reset_index(name='total_count')
    totals_df.to_csv(args.totals_out, sep='\t', index=False)
    print(f"✂︎ Wrote combo totals ➜ {args.totals_out}")

    # 2) Per‐chromosome summary
    combo_order = list(totals.index)
    df['combo'] = pd.Categorical(df['combo'],
                                 categories=combo_order,
                                 ordered=True)
    df = df.sort_values(['combo','chromosome'])
    df.to_csv(args.summary_out, sep='\t', index=False)
    print(f"✂︎ Wrote per‐chromosome summary ➜ {args.summary_out}")

    # 3) Heatmap
    heat = df.pivot(index='combo',
                    columns='chromosome',
                    values='count').fillna(0)
    heat = heat.reindex(combo_order)

    plt.figure(figsize=(8, max(2, 0.4*len(combo_order))))
    im = plt.imshow(heat, aspect='auto')
    plt.colorbar(im, label='Fragment count')
    plt.yticks(range(len(heat.index)), heat.index)
    plt.xticks(range(len(heat.columns)), heat.columns, rotation=90)
    plt.xlabel("Chromosome"); plt.ylabel("Enzyme combo")
    plt.title(f"ddRAD fragments {args.min}–{args.max} bp")
    plt.tight_layout()
    plt.savefig(args.heatmap_out, dpi=150)
    print(f"✂︎ Wrote heatmap ➜ {args.heatmap_out}")

if __name__ == "__main__":
    main()
