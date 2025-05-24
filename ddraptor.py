#!/usr/bin/env python3
import sys
import re
import argparse
import multiprocessing as mp
from functools import partial
from itertools import product
from typing import List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData
import pyahocorasick
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

def expand_iupac(motif: str) -> List[str]:
    """
    Expand an IUPAC‐encoded DNA motif into all possible concrete sequences.
    """
    table = IUPACData.ambiguous_dna_values
    pools = [table.get(b, b) for b in motif]
    return [''.join(p) for p in product(*pools)]


def parse_enzymes(path: str):
    """
    Read an enzyme TSV: each non‐blank, non‐# line is
      enzyme_name <whitespace> motif_with_caret
    Returns list of (name, motif, cut_index).
    """
    enzymes = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            parts = re.split(r'\s+', ln)
            if len(parts) < 2 or '^' not in parts[1]:
                print(f"Warning: skipping malformed enzyme line: {ln}", file=sys.stderr)
                continue
            name, raw = parts[0], parts[1]
            cut_idx = raw.index('^')
            motif   = raw.replace('^', '')
            enzymes.append((name, motif, cut_idx))
    return enzymes


def parse_combos(path: str):
    """
    Read a combo TSV: each non‐blank, non‐# line is
      combo_name <whitespace> enzymeA,enzymeB
    Returns dict combo_name -> (enzymeA, enzymeB).
    """
    combos = {}
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            parts = re.split(r'\s+', ln)
            if len(parts) < 2 or ',' not in parts[1]:
                print(f"Warning: skipping malformed combo line: {ln}", file=sys.stderr)
                continue
            combo, pair = parts[0], parts[1]
            a, b = pair.split(',', 1)
            combos[combo] = (a, b)
    return combos


def build_automaton(enzymes):
    """
    Build a single Aho–Corasick automaton over all
    IUPAC‐expanded motifs and their reverse complements.
    Payload per pattern: (enzyme_idx, cut_idx, pattern_length).
    """
    A = pyahocorasick.Automaton()
    for idx, (name, motif_iupac, cut_idx) in enumerate(enzymes):
        for pat in expand_iupac(motif_iupac):
            L = len(pat)
            # forward
            A.add_word(pat, (idx, cut_idx, L))
            # reverse complement
            rc = str(Seq(pat).reverse_complement())
            if rc != pat:
                rc_cut = L - cut_idx
                A.add_word(rc, (idx, rc_cut, L))
    A.make_automaton()
    return A


def count_ddrad_fragments_between(
    cuts1: List[int],
    cuts2: List[int],
    min_len: int,
    max_len: int
) -> int:
    """
    Count ddRAD fragments flanked by cuts in cuts1 and cuts2
    using the original “start→end, end→next start” logic,
    but compute counts directly.
    """
    if not cuts1 or not cuts2:
        return 0

    # Determine which enzyme cuts come first
    A, B = (cuts1, cuts2) if cuts1[0] < cuts2[0] else (cuts2, cuts1)
    i = j = count = 0
    while i < len(A) and j < len(B):
        # advance A until A[i] > B[j]
        while i < len(A) and A[i] <= B[j]:
            i += 1
        if i >= len(A):
            break
        # fragment A[i-1] -> B[j]
        length = B[j] - A[i-1]
        if min_len <= length <= max_len:
            count += 1

        # advance B until B[j] > A[i]
        while j < len(B) and B[j] <= A[i]:
            j += 1
        if j >= len(B):
            break
        # fragment B[j-1] -> A[i]
        length = A[i] - B[j-1]
        if min_len <= length <= max_len:
            count += 1

    return count


def count_ddrad_fragments(record, enzymes, automaton, combos, min_len, max_len):
    """
    For a single SeqRecord:
    - find all cut sites,
    - collect per-enzyme cut lists,
    - count ddRAD fragments for each combo.
    Returns list of (combo, chromosome, count).
    """
    seq = str(record.seq)
    # collect cuts per enzyme
    cuts = { name: [] for name, _, _ in enzymes }
    for end_pos, (idx, cut_idx, L) in automaton.iter(seq):
        start_pos = end_pos - L + 1
        pos = start_pos + cut_idx
        cuts[enzymes[idx][0]].append(pos)
    # sort each enzyme’s cut list
    for lst in cuts.values():
        lst.sort()

    results = []
    for combo, (e1, e2) in combos.items():
        cnt = count_ddrad_fragments_between(cuts[e1], cuts[e2], min_len, max_len)
        results.append((combo, record.id, cnt))
    return results


def main():
    p = argparse.ArgumentParser(
        description="ddRAD simulator with fast pairing logic"
    )
    p.add_argument('enzymes_tsv', help="TSV: enzyme_name <tab> motif^cut")
    p.add_argument('combos_tsv',   help="TSV: combo_name <tab> enzA,enzB")
    p.add_argument('fasta',        help="Reference FASTA")
    p.add_argument('--min',  type=int, required=True, help="Min fragment length")
    p.add_argument('--max',  type=int, required=True, help="Max fragment length")
    p.add_argument('--processes', type=int, default=mp.cpu_count(),
                   help="Number of worker processes")
    p.add_argument('--totals-out',  default='ddrad_totals.tsv',
                   help="Combo totals TSV path")
    p.add_argument('--summary-out', default='ddrad_summary.tsv',
                   help="Per-chromosome summary TSV path")
    p.add_argument('--heatmap-out', default='ddrad_heatmap.png',
                   help="Heatmap PNG path")
    args = p.parse_args()

    # load
    enzymes = parse_enzymes(args.enzymes_tsv)
    combos  = parse_combos(args.combos_tsv)
    autom   = build_automaton(enzymes)

    # worker
    worker = partial(
        count_ddrad_fragments,
        enzymes=enzymes,
        automaton=autom,
        combos=combos,
        min_len=args.min,
        max_len=args.max
    )

    # parallel processing
    records = list(SeqIO.parse(args.fasta, "fasta"))
    pool    = mp.Pool(args.processes)
    all_rows = []
    for rec_res in tqdm(pool.imap(worker, records),
                        total=len(records), desc="Digesting"):
        all_rows.extend(rec_res)
    pool.close(); pool.join()

    # build DataFrame
    df = pd.DataFrame(all_rows, columns=['combo','chromosome','count'])

    # combo totals
    totals = df.groupby('combo')['count'].sum().sort_values(ascending=False)
    totals_df = totals.reset_index(name='total_count')
    totals_df.to_csv(args.totals_out, sep='\t', index=False)
    print(f"✂︎ Wrote combo totals ➜ {args.totals_out}")

    # per‐chromosome summary
    combo_order = list(totals.index)
    df['combo'] = pd.Categorical(df['combo'], categories=combo_order, ordered=True)
    df = df.sort_values(['combo','chromosome'])
    df.to_csv(args.summary_out, sep='\t', index=False)
    print(f"✂︎ Wrote per-chromosome summary ➜ {args.summary_out}")

    # heatmap
    heat = df.pivot(index='combo', columns='chromosome', values='count').fillna(0)
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
