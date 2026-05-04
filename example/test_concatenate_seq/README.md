# Test Examples for concatenate_seq.py

This directory contains test data and examples for testing the `concatenate_seq.py` script functionality. The script concatenates multiple NEXUS alignment files and produces a single FASTA file with sequences ordered alphabetically by taxa names.

## Directory Structure

```
test_concatenate_seq/
├── README.md           # This file
├── rep1/               # Test replicate 1
│   ├── gene1.nex      # Gene 1 alignment (7 taxa)
│   ├── gene2.nex      # Gene 2 alignment (6 taxa) 
│   └── gene3.nex      # Gene 3 alignment (4 taxa)
└── rep2/              # Test replicate 2
    ├── gene1.nex      # Gene 1 alignment (6 taxa)
    ├── gene2.nex      # Gene 2 alignment (4 taxa)
    └── gene3.nex      # Gene 3 alignment (4 taxa)
```

## Test Data Features

The test data is designed to test several important features:

1. **Alphabetical Sorting**: Taxa are named with patterns like `A_0`, `A_1`, `A_2`, `B_0`, `C_0`, `D_0`, `D_1` to test proper alphabetical sorting by letter first, then by number.

2. **Missing Taxa Handling**: Different gene files contain different sets of taxa to test the script's ability to handle missing taxa by filling gaps with `-` characters.

3. **Multiple Replicates**: Two replicate directories test batch processing capabilities.

## Sample Taxa in Test Data

- **rep1**: Contains taxa `A_0`, `A_1`, `A_2`, `B_0`, `C_0`, `C_1`, `D_0`, `D_1` (not all present in every gene)
- **rep2**: Contains taxa `A_0`, `A_1`, `B_0`, `C_0`, `D_0` (not all present in every gene)

## How to Run the Tests
### Basic Command Line Usage

Navigate to the simulation-reptiles root directory and run:

```bash
# Test replicate 1
python scripts/concatenate_seq.py \
    example/test_concatenate_seq/rep1 \
    example/test_concatenate_seq/ \
    1

# Test replicate 2  
python scripts/concatenate_seq.py \
    example/test_concatenate_seq/rep2 \
    example/test_concatenate_seq/ \
    2
```
These example commands demonstrate how to run the `concatenate_seq.py` script for different test replicates.  
- The first argument specifies the input directory for each replicate (`rep1` or `rep2`).  
- The second argument is the output directory where results will be saved.  
- The third argument indicates the replicate number.

## Expected Output

After running the commands, you should see these output files:

```
test_concatenate_seq/
├── concate_alignment_rep1.fasta    # Concatenated sequences from rep1
└── concate_alignment_rep2.fasta    # Concatenated sequences from rep2
```

### Sample Output Format

Each FASTA file will contain sequences ordered alphabetically:

```
>A_0
ACTCTCTCGACCTCTCTCGAACTCTCTCGA
>A_1
CCTCTCTCGACCTCTCTCGCACTCTCTCCC
>A_2
CCTCTCTCCACCTCTCTCGT----------
>B_0
ACTCTCTCGACCTCTCTCGAACTCTCTCGA
>C_0
ACTCTCTCGGCCTTTCTCGGACTCTCTCGG
>C_1
------------------------ACTCTCTCGT
>D_0
ACTCTCTTTTCCTCTCTTTT----------
>D_1
ACTCTCTTTG------------------------
```

**Note:** Sequences with `-` characters indicate missing taxa in some gene files, which the script handles by padding with gaps.