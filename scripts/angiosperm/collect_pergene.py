#!/usr/bin/env python3
"""
Collect per-gene IQ-TREE runs into the files used by ASTRAL and by
scripts/speciestree_basal_angio.jl. Replaces `separate-boot-bygene.jl`
from the reptiles repo.

Why per-gene runs and not `iqtree -S alignments` (one partitioned run):
with -S, when a gene lacks some taxa, IQ-TREE 3.1.3 writes that gene's
bootstrap trees (-wbtl, loci.ufboot) with wrong tip names (a missing
taxon's name replaces a present one), which corrupts the bootstrap support
computed by ASTRAL. The ML trees are correct. Running IQ-TREE once per
gene avoids the problem:
  cd data/a353/<name>
  for f in alignments/*.fasta; do g=$(basename $f .fasta)
    iqtree -s $f -st DNA -B 1000 -wbtl -m MFP -mset HKY -mrate G -mfreq F \
           -seed 1 -T 1 --prefix iqtree_pergene/$g
  done
Then (from data/a353/<name>):  python3 collect_pergene.py
which writes
  loci.txt              gene names, in the order used below
  iqtree/loci.treefile  ML gene trees, one per line, same order
  iqtree/bootstrap/<gene>.ufboot   bootstrap trees per gene
  BSlistfiles           one bootstrap-file path per line, same order
  seqgen_params.csv     per-gene HKY+F+G parameters (kappa, alpha, freqs)
                        from the .iqtree reports, to be fitted as in
                        notes/choice-seqgen-parameters.md
Paths in BSlistfiles are relative to the current folder: run ASTRAL from
the same folder:
  astral -i iqtree/loci.treefile -b BSlistfiles -r 1000 -o astral/species.tre
"""
import glob
import os
import re

src = "iqtree_pergene"
genes = sorted(os.path.basename(f)[:-9] for f in glob.glob(f"{src}/*.treefile"))
os.makedirs("iqtree/bootstrap", exist_ok=True)
os.makedirs("astral", exist_ok=True)

with open("loci.txt", "w") as lo, open("iqtree/loci.treefile", "w") as tf, \
        open("BSlistfiles", "w") as bl, open("seqgen_params.csv", "w") as sp:
    sp.write("gene,kappa,alpha,fA,fC,fG,fT\n")
    for g in genes:
        lo.write(g + "\n")
        tf.write(open(f"{src}/{g}.treefile").read().strip() + "\n")
        boot = f"iqtree/bootstrap/{g}.ufboot"
        with open(boot, "w") as o:
            o.write(open(f"{src}/{g}.ufboot").read())
        bl.write(boot + "\n")
        rep = open(f"{src}/{g}.iqtree").read()
        kappa = re.search(r"A-G:\s*([0-9.]+)", rep)
        alpha = re.search(r"Gamma shape alpha:\s*([0-9.]+)", rep)
        fr = re.findall(r"pi\([ACGT]\)\s*=\s*([0-9.]+)", rep)
        if kappa and alpha and len(fr) == 4:
            sp.write(f"{g},{kappa.group(1)},{alpha.group(1)},"
                     + ",".join(fr) + "\n")
        else:
            print(f"warning: no HKY+F+G parameters parsed for {g}")
print(f"collected {len(genes)} genes")
