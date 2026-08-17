"""
Convert a concatenated FASTA alignment to sequential PHYLIP for PhyNEST.

Uses Biopython (SeqIO/AlignIO), the same library used by
concatenate_seq.py. The "phylip-sequential" writer keeps each full
sequence on a single line, which PhyNEST's readPhylip requires
(wrapped/interleaved PHYLIP is not supported by PhyNEST).

Tip names like "A_0" are simplified to species letters ("A").
PhyNEST has no individual-to-species mapping for its network search,
so when --n_inds >= 2 only the first individual ("_0") of each
species is kept (see phynest_1rep.jl).

### Example Usage:
python scripts/fasta2phylip.py \
    -i rep01/seqgenfolder/concate_alignment_rep01.fasta \
    -o rep01/phynestfolder/rep01.phy --n_inds 1
"""

import argparse
import sys

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment


def convert(fasta, phylip, n_inds):
    """Read FASTA, subsample/rename tips, write sequential PHYLIP."""
    records = list(SeqIO.parse(fasta, "fasta"))
    if n_inds >= 2:
        records = [r for r in records if r.id.endswith("_0")]
    for r in records:
        r.id = r.id.split("_")[0]
        r.description = ""

    names = [r.id for r in records]
    if len(set(names)) != len(names):
        sys.exit(f"ERROR: duplicate names after simplifying: {names}")
    if any(len(n) > 10 for n in names):
        # phylip-sequential truncates names longer than 10 characters
        sys.exit(f"ERROR: name longer than 10 characters in: {names}")
    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        sys.exit(f"ERROR: unequal sequence lengths: {lengths}")

    alignment = MultipleSeqAlignment(records)
    AlignIO.write(alignment, phylip, "phylip-sequential")
    print(f"Wrote {len(names)} taxa x {lengths.pop()} bp to {phylip}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="FASTA to sequential PHYLIP for PhyNEST")
    parser.add_argument("-i", "--input", required=True,
                        help="input FASTA alignment")
    parser.add_argument("-o", "--output", required=True,
                        help="output PHYLIP file")
    parser.add_argument("--n_inds", type=int, default=1,
                        help="individuals per species; if >= 2, keep "
                             "only the first ('_0') individual")
    args = parser.parse_args()
    convert(args.input, args.output, args.n_inds)
