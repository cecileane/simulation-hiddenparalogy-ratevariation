#!/usr/bin/env python3
"""
Subset the 1KP capstone FNA2AA single-copy gene alignments to a chosen set of
8 samples (4-letter 1KP codes), keeping only genes present in all 8 taxa and
removing codon columns that are gaps in all 8 taxa.

Input (downloaded from CyVerse, see notes/angiosperm_speciestree_plan.md):
  data/onekp/FNA2AA/genes/<gene>/FNA2AA-upp-masked.fasta   (374 genes)
  data/onekp/annotations.csv                                (code -> species)

Usage (from repo root):
  python3 scripts/angiosperm/subset_onekp.py --name ginkgo \
      --codes SGTW URDJ PZRT WTKZ VZCI WZFE XQWC NPND
Writes data/onekp/subset_<name>/alignments/<gene>.fasta and taxa.csv
"""
import argparse, csv, glob, os, statistics

p = argparse.ArgumentParser()
p.add_argument("--name", required=True)
p.add_argument("--codes", nargs=8, required=True, help="8 codes, in order A..H")
p.add_argument("--onekp", default="data/onekp")
a = p.parse_args()

ann = {r["Code"]: r["Species"] for r in csv.DictReader(open(f"{a.onekp}/annotations.csv"))}
files = sorted(f for f in glob.glob(f"{a.onekp}/FNA2AA/genes/*/*.fasta") if os.path.isfile(f))

def readfa(f):
    d, k = {}, None
    for l in open(f):
        l = l.strip()
        if l.startswith(">"):
            k = l[1:]; d[k] = []
        else:
            d[k].append(l)
    return {k: "".join(v) for k, v in d.items()}

out = f"{a.onekp}/subset_{a.name}/alignments"
os.makedirs(out, exist_ok=True)
lens = []
for f in files:
    seqs = readfa(f)
    if not all(c in seqs for c in a.codes):
        continue
    gene = os.path.basename(os.path.dirname(f))
    L = len(seqs[a.codes[0]])
    assert L % 3 == 0, f"{gene}: alignment length not a multiple of 3"
    keep = [j for j in range(0, L, 3) if any(seqs[c][j:j+3].replace("-", "") for c in a.codes)]
    with open(f"{out}/{gene}.fasta", "w") as o:
        for c in a.codes:
            o.write(f">{c}\n" + "".join(seqs[c][j:j+3] for j in keep) + "\n")
    lens.append(3 * len(keep))

with open(f"{a.onekp}/subset_{a.name}/taxa.csv", "w") as o:
    o.write("letter,code,species\n")
    for L, c in zip("ABCDEFGH", a.codes):
        o.write(f"{L},{c},{ann.get(c, 'NA')}\n")
print(f"{a.name}: {len(lens)} genes with all 8 taxa; length min/median/max = "
      f"{min(lens)}/{statistics.median(lens)}/{max(lens)}; total {sum(lens)} bp")
