# ddRaptor

**ddRaptor** is a ddRAD fragment simulator that uses Aho–Corasick pattern matching to find cut sites for enzyme pairs, counts the fragments within a size range, and generates summary tables and a heatmap of fragment counts per chromosome.

---

## Dependencies

* Python 3.6+
* [Biopython](https://biopython.org/)
* [pyahocorasick](https://pypi.org/project/pyahocorasick/)
* [tqdm](https://pypi.org/project/tqdm/)
* [pandas](https://pypi.org/project/pandas/)
* [matplotlib](https://matplotlib.org/)

Install them via:

```bash
pip install biopython pyahocorasick tqdm pandas matplotlib
```

---

## Installation

1. Clone or download this repository.

2. Make sure you have the dependencies installed (see above).

3. Ensure the script is executable:

   ```bash
   chmod +x ddraptor.py
   ```

4. (Optional) Copy into your `PATH`:

   ```bash
   cp ddraptor.py /usr/local/bin/ddraptor
   ```

---

## Input File Formats

### 1. `enzymes.tsv`

Tab-delimited list of enzymes and their IUPAC motifs, with a caret (`^`) marking the cut position:

```tsv
# name      motif_with_caret
EcoRI       G^AATTC
BamHI       G^GATCC
HindIII     A^AGCTT
SphI        GCATG^C
PstI        CTGCA^G
CviQI       G^TAC
NsiI        ATGCA^T
```

### 2. `combos.tsv`

Tab-delimited list of ddRAD enzyme pairs (combos):

```tsv
# combo_name    enzymeA,enzymeB
Combo1         EcoRI,BamHI
Combo2         EcoRI,HindIII
Combo3         BamHI,HindIII
Combo4         SphI,PstI
Combo5         EcoRI,PstI
Combo6         CviQI,NsiI
```

### 3. Reference FASTA

Any multi-FASTA file with chromosome or contig sequences:

```fasta
>chr1
ACGTACGT...
>chr2
TTGACGTA...
...
```

---

## Usage

```bash
python ddraptor.py \
  <enzymes.tsv> <combos.tsv> <reference.fasta> \
  --min <MIN_LENGTH> --max <MAX_LENGTH> [options]
```

### Required arguments

* `enzymes.tsv` – Path to enzyme definitions.
* `combos.tsv` – Path to combo definitions.
* `reference.fasta` – Path to your multi-FASTA reference.
* `--min`  Minimum fragment length (inclusive).
* `--max`  Maximum fragment length (inclusive).

### Optional arguments

* `--processes`   Number of parallel worker processes (default: number of CPU cores).
* `--totals-out`   Path for combo totals TSV (default: `ddrad_totals.tsv`).
* `--summary-out`  Path for per-chromosome summary TSV (default: `ddrad_summary.tsv`).
* `--heatmap-out`  Path for heatmap PNG (default: `ddrad_heatmap.png`).

---

## Output

1. **Combo totals TSV** (`--totals-out`)
   Columns: `combo` | `total_count`

   ```tsv
   combo    total_count
   Combo3   12345
   Combo1    9876
   ...
   ```

2. **Per-chromosome summary TSV** (`--summary-out`)
   Columns: `combo` | `chromosome` | `count`

   ```tsv
   combo    chromosome    count
   Combo3   chr1          2345
   Combo3   chr2          1987
   ...
   Combo1   chr1          1234
   ...
   ```

3. **Heatmap PNG** (`--heatmap-out`)
   A matrix plot of fragment counts (`count`) with:

   * **Rows** = enzyme combinations, sorted by `total_count` descending
   * **Columns** = chromosomes/contigs in the FASTA

---

## Examples

### Default outputs

```bash
python ddraptor.py enzymes.tsv combos.tsv genome.fasta \
  --min 200 --max 600
```

Produces:

* `ddrad_totals.tsv`
* `ddrad_summary.tsv`
* `ddrad_heatmap.png`

### Custom output paths

```bash
python ddraptor.py enzymes.tsv combos.tsv genome.fasta \
  --min 150 --max 600 --processes 8 \
  --totals-out=my_totals.tsv \
  --summary-out=my_summary.tsv \
  --heatmap-out=my_heatmap.png
```

---

## Tips & Notes

* **IUPAC & reverse strands:**
  The script auto-expands ambiguous IUPAC motifs and searches both forward and reverse complements.
* **Performance:**
  Uses Aho–Corasick plus multiprocessing; scales linearly with CPU cores and number of contigs.
* **Error handling:**
  Malformed lines in TSVs are skipped with a warning to stderr.

---

> **License:** MIT \
> **Author:** Georgios Kousis Tsampazis \
> **Contact:** [georgekousis6@gmail.com](mailto:georgekousis6@gmail.com)
