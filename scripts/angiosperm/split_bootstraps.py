#!/usr/bin/env python3
"""
Mimics `separate-boot-bygene.jl` from the reptiles repo.

After `iqtree -S alignments --prefix iqtree/loci -B 1000 -wbtl`, IQ-TREE writes
  iqtree/loci.treefile : one ML tree per partition (gene), in partition order
  iqtree/loci.ufboot   : all bootstrap trees, 1000 per partition, same order
This script splits loci.ufboot into one file per gene and writes BSlistfiles
(one bootstrap-file path per line, same order as loci.treefile) for
  astral -i iqtree/loci.treefile -b BSlistfiles -r 1000

Usage (from the subset folder):  python3 split_bootstraps.py iqtree 1000
"""
import os, re, sys

prefix_dir = sys.argv[1] if len(sys.argv) > 1 else "iqtree"
B = int(sys.argv[2]) if len(sys.argv) > 2 else 1000

# partition order from the best_scheme.nex file (charset lines)
names = []
for l in open(f"{prefix_dir}/loci.best_scheme.nex"):
    m = re.match(r"\s*charset\s+(\S+)\s*=", l)
    if m:
        names.append(m.group(1))
ml = [l for l in open(f"{prefix_dir}/loci.treefile") if l.strip()]
boot = [l for l in open(f"{prefix_dir}/loci.ufboot") if l.strip()]
assert len(ml) == len(names), f"{len(ml)} ML trees vs {len(names)} partitions"
assert len(boot) == B * len(names), f"{len(boot)} bootstrap trees, expected {B*len(names)}"

os.makedirs(f"{prefix_dir}/bootstrap", exist_ok=True)
with open("BSlistfiles", "w") as lst:
    for i, n in enumerate(names):
        fn = f"{prefix_dir}/bootstrap/{n}.ufboot"
        with open(fn, "w") as o:
            o.writelines(boot[i*B:(i+1)*B])
        lst.write(fn + "\n")
print(f"wrote {len(names)} bootstrap files and BSlistfiles")
