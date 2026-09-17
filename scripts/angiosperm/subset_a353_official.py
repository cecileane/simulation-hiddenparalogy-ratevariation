#!/usr/bin/env python3
"""
Build an 8-taxon per-gene dataset from the official per-gene alignments of
the Kew Tree of Life Explorer release (Angiosperms353; Zuntini et al. 2024,
Nature 629:843), so that no alignment step is done by us.

Source (anonymous):
  https://sftp.kew.org/pub/paftol/current_release/fasta/alignments/
  <gene>.dna.aln.fasta : 353 files, 20-27 MB each (~20,000 samples per
  gene), with headers like
  >4527 Gene_Name:PUR2 Species:Amborella_trichopoda Repository:INSDC
        Sequence_ID:GCF_000471905.2
Each file is downloaded, the 8 sequences are pulled out by Sequence_ID,
columns that are gaps in all 8 are removed, and the download is deleted,
so only ~2 MB of extracted alignments remain on disk.

Genes are kept if at least --min_taxa of the 8 samples are present
(default 4; IQ-TREE and ASTRAL handle missing taxa, and the reptile UCE
loci did not all contain all taxa either). Presence/absence per gene is
written to gene_occupancy.csv. The run is resumable: genes already listed
in gene_occupancy.csv are skipped.

Usage (from repo root):
  python3 scripts/angiosperm/subset_a353_official.py --name pinus \\
      --samples CODE=Sequence_ID ... (8 entries, in order A..H, A = outgroup)
Writes data/a353/<name>/alignments/<gene>.fasta, taxa.csv, gene_occupancy.csv
"""
import argparse
import os
import re
import statistics
import subprocess

BASE = "https://sftp.kew.org/pub/paftol/current_release"
p = argparse.ArgumentParser()
p.add_argument("--name", required=True)
p.add_argument("--samples", nargs=8, required=True,
               help="CODE=Sequence_ID, 8 entries, A..H")
p.add_argument("--tmp", default="/tmp")
p.add_argument("--min_taxa", type=int, default=4)
a = p.parse_args()

codes, ids = [], []
for s in a.samples:
    c, i = s.split("=", 1)
    codes.append(c)
    ids.append(i)

# gene list from the release manifest
man = subprocess.run(["curl", "-sS", f"{BASE}/gene_manifest.txt"],
                     capture_output=True, text=True, check=True).stdout
genes = [l.split("\t")[0] for l in man.splitlines()[1:] if l.strip()]
out = f"data/a353/{a.name}"
os.makedirs(f"{out}/alignments", exist_ok=True)

# resumable: genes already listed in gene_occupancy.csv are skipped
occfile = f"{out}/gene_occupancy.csv"
done = set()
if os.path.exists(occfile):
    done = {l.split(",")[0] for l in open(occfile).read().splitlines()[1:]}
    occ = open(occfile, "a")
else:
    occ = open(occfile, "w")
    occ.write("gene," + ",".join(codes) + "\n")

species = {}
lens = []
for k, g in enumerate(genes, 1):
    if g in done:
        continue
    tmp = os.path.join(a.tmp, f"{g}.dna.aln.fasta")
    subprocess.run(["curl", "-sS", "--retry", "5", "--retry-delay", "10",
                    "--retry-all-errors", "-o", tmp,
                    f"{BASE}/fasta/alignments/{g}.dna.aln.fasta"],
                   check=True)
    seqs, cur = {}, None
    with open(tmp) as fh:
        for l in fh:
            if l.startswith(">"):
                m = re.search(r"Sequence_ID:(\S+)", l)
                sid = m.group(1) if m else None
                cur = sid if sid in ids else None
                if cur:
                    seqs[cur] = []
                    sp = re.search(r"Species:(\S+)", l)
                    species[cur] = (sp.group(1).replace("_", " ")
                                    if sp else "NA")
            elif cur:
                seqs[cur].append(l.strip())
    os.remove(tmp)
    seqs = {k2: "".join(v) for k2, v in seqs.items()}
    occ.write(g + "," + ",".join("1" if i in seqs else "0" for i in ids)
              + "\n")
    occ.flush()
    present = [i for i in ids if i in seqs]
    if len(present) < a.min_taxa:
        print(f"[{k}/{len(genes)}] {g}: only {len(present)} taxa, skipped",
              flush=True)
        continue
    L = len(seqs[present[0]])
    keep = [j for j in range(L) if any(seqs[i][j] != "-" for i in present)]
    with open(f"{out}/alignments/{g}.fasta", "w") as o:
        for c, i in zip(codes, ids):
            if i in seqs:
                o.write(f">{c}\n" + "".join(seqs[i][j] for j in keep) + "\n")
    lens.append(len(keep))
    miss = [c for c, i in zip(codes, ids) if i not in seqs]
    print(f"[{k}/{len(genes)}] {g}: {len(keep)} columns, {len(present)} taxa"
          + (f", missing {miss}" if miss else ""), flush=True)
occ.close()

if species:  # taxa.csv is written when at least one gene was processed
    with open(f"{out}/taxa.csv", "w") as o:
        o.write("letter,code,species\n")
        for Lt, c, i in zip("ABCDEFGH", codes, ids):
            o.write(f"{Lt},{c},{species.get(i, 'NA')}\n")
n = len(os.listdir(f"{out}/alignments"))
print(f"\n{a.name}: {n} genes with >= {a.min_taxa} samples")
if lens:
    print(f"this run: aligned length min/median/max = {min(lens)}/"
          f"{statistics.median(lens)}/{max(lens)}; total {sum(lens)} bp")
