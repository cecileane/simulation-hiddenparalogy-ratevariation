# Scripts Folder README
This folder contains the main simulation scripts and utilities for the simulation-reptiles project. The workflow consists of three main stages: data simulation, network inference, and graph analysis. 

## Complete Workflow Overview

```
Stage 1                 Stage 2 
simulation_iqtree.jl  → findgraphs.jl
                      → snaq_submit.jl
```

1. **Data Simulation** (`simulation_iqtree.jl`): Generates simulated gene trees handling hidden paralogy, substitution rate variation across genes and/or lineages. `simulation_iqtree.jl` is the major pipeline script which calls the following helper and wrapper scripts:

   - **[`utilities.jl`](utilities.jl)**: Provides functions for tree processing (e.g., filtering gene trees, renaming tips), seed generation, file management, and parameter handling.
   - **[`seq-gen.sh`](seq-gen.sh)**: Bash script that wraps Seq-Gen calls for simulating molecular sequences along gene trees, handling input/output organization.
   - **[`concatenate_seq.py`](concatenate_seq.py)**: Python script for concatenating individual gene alignments into a single multi-gene alignment in FASTA or NEXUS format, used for downstream analysis.
   - **[`iqtree.pl`](iqtree.pl)**: Perl script that automates batch execution of IQ-TRE E and generates species tree using ASTER (astral) across multiple gene alignments, manages parallelization, and collects output files.

The output of data of `simulation_iqtree.jl` is used for the following processes respectively: 
2.1. **Admixture Graph Analysis** (`findgraphs.jl`): Constructs population graphs using AdmixTools. 
2.2. **Network Inference** (`snaq_submit.jl`): Infers phylogenetic networks using SNaQ. 


## Stage 1: Data Simulation (simulation_iqtree.jl)

### Overview
[`simulation_iqtree.jl`](simulation_iqtree.jl) implements a pipeline that:

1. _Simulates gene trees_ using SimPhy with hidden paralogy modeling
2. _Generates molecular sequences_ using Seq-Gen 
3. _Estimates gene trees_ using IQ-TREE
4. _Reconstructs species trees_ using ASTRAL

### Software Dependencies

The pipeline calls the following external software. Installation instructions can be found in [software_installation.sh](../notebook/software_installation.sh).

- _[SimPhy](https://github.com/adamallo/SimPhy)_: Simulates gene trees along a species tree with gene duplication/loss
- _[Seq-Gen](https://github.com/rambaut/Seq-Gen)_: Simulates molecular sequence evolution along gene trees
- _[IQ-TREE](http://www.iqtree.org/)_: Maximum likelihood phylogenetic inference
- _[ASTRAL](https://github.com/smirarab/ASTRAL)_: Species tree estimation from gene trees

All executables should be placed in the [`executables/`](../executables/) folder.

### Detailed Workflow

#### 1. Gene Tree Simulation (SimPhy + Hidden Paralogy Modeling)

**Purpose**: Simulate gene trees that may contain hidden paralogy events.

**Process**:
- Uses SimPhy to generate gene trees along a fixed species tree topology
- Supports gene duplication (`--dup_rate`) and loss (`--loss_rate`) rates
- Implements rate variation across genes (`G`), lineages (`L`), or both (`GL`) or no variation (`N`). 
- **Hidden Paralogy Logic**: After SimPhy generates trees, the pipeline modifies them to simulate hidden paralogy:
  
  **Trees are EXCLUDED if they have**:
  - More than one gene copy per individual (obvious paralogy, not hidden)
  - ≤3 taxa remaining after gene loss (phylogenetically uninformative)
  
  **Trees are RETAINED if they have**:
  - 0-1 gene copy per individual and ≥4 taxa (potential hidden paralogy)
  - Gene copy IDs are removed from tip names, keeping only species + individual IDs

**Iterative Re-running**: 
- If insufficient trees are generated, SimPhy is re-run up to `max_iteration_simphy` times
- Pipeline continues if `≥ min_gene_proportion * n_genes` trees are obtained
- Warns if final count is below minimum threshold

#### 2. Molecular Sequence Simulation (Seq-Gen)

**Purpose**: Generate DNA alignments along the filtered gene trees.

**Process**:
- Uses each retained gene tree as input to Seq-Gen
- Simulates 1000bp sequences per gene (based on empirical UCE data)
- Uses HKY substitution model with parameters from real reptile data:
  - Transition/transversion ratio (kappa) = 4.143
  - Base frequencies: A=0.316, T=0.182, G=0.183, C=0.319
  - Gamma-distributed rate variation across sites

#### 3. Gene Tree Estimation (IQ-TREE)

**Purpose**: Reconstruct gene trees from simulated sequences to test inference accuracy.

**Process**:
- Runs IQ-TREE on each simulated alignment
- Uses HKY+G model (matching simulation parameters)
- No bootstrap analysis (for speed)
- Processes multiple genes per replicate in parallel
- Uses [`iqtree.pl`](iqtree.pl) wrapper script for batch processing

#### 4. Species Tree Reconstruction (ASTRAL)

**Purpose**: Infer species trees from estimated gene trees using coalescent methods.

**Process**:
- Concatenates all gene trees from IQ-TREE into a single file
- Creates mapping file for multiple individuals per species (if `n_inds > 1`)
- Runs ASTRAL to estimate species tree with local posterior probabilities
- Uses wASTRAL (weighted ASTRAL) implementation

### Key Parameters

#### Rate Variation Options (`--ratevar`):
- `N`: No rate variation
- `G`: Gene-specific rate heterogeneity (log-normal distribution)
- `L`: Lineage-specific rate variation (branch-specific multipliers)  
- `GL` or `G*L`: Both gene and lineage rate variation

#### Seed Management

- **Deterministic seeding system**: A master seed is generated for each unique parameter combination (e.g., duplication/loss rates, rate variation, number of individuals, etc) by assigning each parameter a unique prime number and multiplying the prime numbers together.
- **Stage-specific seeds**: Each major software component (SimPhy, Seq-Gen, IQ-TREE, ASTRAL) receives its own seed, deterministically derived from the master seed and the replicate number.
- **Replicate-level reproducibility**: For each replicate, seeds are generated and stored in a text file (e.g., `random_seed_simphy.txt`, `random_seed_findgraphs.txt`, etc) within the replicate folder. This file records the seeds used for all tools in that replicate.
- **Downstream consistency**: The same seed derivation logic is used in downstream scripts (`snaq_submit.jl`, `findgraphs.jl`) to ensure reproducibility across all pipeline stages.
- **Implementation**: See [`utilities.jl`](utilities.jl) for the `generate_seeds()` and related functions that handle seed calculation, assignment, and file output.

### Output Structure

For each parameter combination, creates folder: `output/DUP{rate}-LOS{rate}-RV{var}-N_ind{n}/`

Within each replicate folder (`rep001/`, `rep002/`, etc.):
```
rep001/
├── genetrees_simphy/     # Raw SimPhy output + configuration files
├── genetrees_singlecopy/ # Filtered trees after hidden paralogy  
├── seqgenfolder/         # Seq-Gen alignments (nexus + concatenated fasta)
├── iqtreefolder/         # IQ-TREE gene trees and logs
└── astralfolder/         # ASTRAL species tree
```

### Parallel Processing

The pipeline uses Julia's `Distributed` package for parallelization:
- Each replicate processed independently
- Launch with: `julia -p N simulation_iqtree.jl --args...` (where N = number of cores)

### Logging and Tracking

**Simulation tracking** (`Simphy_{params}.csv`):
- Records number of gene trees generated per replicate
- Tracks trees removed due to repeated/insufficient taxa
- Number of SimPhy iterations required

**Warning logs** (`{params}.log`):
- Reports replicates with insufficient gene trees
- Used by downstream analysis scripts ([`snaq_submit.jl`](snaq_submit.jl), [`findgraphs.jl`](findgraphs.jl))

**Timing information** (`arguments-{params}.log`):
- Detailed runtime for each pipeline stage
- Complete parameter and seed information
- Concatenation warnings and messages

## Related Scripts

- [`utilities.jl`](utilities.jl): Helper functions for tree processing, seed generation, file management
- [`iqtree.pl`](iqtree.pl): Perl wrapper for batch IQ-TREE processing
- [`seq-gen.sh`](seq-gen.sh): Bash wrapper for Seq-Gen execution
- [`concatenate_seq.py`](concatenate_seq.py): Python script for sequence concatenation
- [`speciestree.jl`](speciestree.jl): Species tree preparation and parameter estimation

## Usage Example

```bash
# Run with 8 cores, simulate hidden paralogy scenario
julia -p 8 scripts/simulation_iqtree.jl \
  --dup_rate 0.01 \
  --loss_rate 0.005 \
  --ratevar GL \
  --n_reps 100 \
  --n_genes 1000 \
  --n_inds 2
```

This creates a simulation with moderate gene duplication/loss rates, both gene and lineage rate variation, 100 replicates of 1000 genes each, with 2 individuals per species.

## Stage 2.1: Admixture Graph Analysis (findgraphs.jl)

### Overview

[`findgraphs.jl`](findgraphs.jl) analyzes the simulated datasets from `simulation_iqtree.jl` using population graph methods to detect admixture patterns and population structure. This script processes the concatenated sequence alignments to construct population graphs and compute admixture statistics using AdmixTools.

### Scripts Called by findgraphs.jl

- [`findgraphs_1rep.jl`](findgraphs_1rep.jl): Core processing script that handles individual replicates
- [`utilities.jl`](utilities.jl): Helper functions for file management and data processing
- Additional utility scripts for data format conversion and statistical analysis

### Software Dependencies

The pipeline calls the following external software:

- **[AdmixTools](https://github.com/DReichLab/AdmixTools)**: Population genetic analysis toolkit
- **[snp-sites](https://github.com/sanger-pathogens/snp-sites)**: Tool for extracting SNP sites from multi-FASTA alignments, used to convert concatenated alignments into SNP matrices suitable for downstream population genetic analysis.

All executables should be placed in the [`executables/`](../executables/) folder.

### Detailed Workflow

#### 1. Input Data Processing

**Purpose**: Convert simulated sequence data into formats suitable for population genetic analysis.

**Process**:
- Reads concatenated FASTA alignments from `seqgenfolder/` and `concatenated_seq.py` (output of `simulation_iqtree.jl`)
- Converts sequences to eigenstrat format (.bed/.bim/.fam files)
- Handles multiple individuals per species when `n_inds > 1`
- Filters sites for minimum allele frequency and missing data thresholds
- Creates population assignment files for downstream analysis

#### 2. Population Graph Construction via findgraphs_1rep.jl

**Purpose**: Process individual replicates to build population graphs and compute admixture statistics.

**How `findgraphs_1rep.jl` is Called**:
```julia
# For each replicate directory
for rep_dir in replicate_directories
   # Call findgraphs_1rep.R with specific parameters
   run(`Rscript findgraphs_1rep.R \
      --input_dir $rep_dir/seqgenfolder \
      --output_dir $rep_dir/admixtools_output \
      --prefix concatenated \
      --num_admix 1 \
      --stop_gen 100 \
      --outgroup A \
      --runs 100 \
      --blgsize 300 \
      --rep_id $rep_id \
      --seed_file_path $seed_file \
      --output_graph_suffix _graphs.rds \
      --output_f2_suffix _f2.rds \
      --output_summary_table_suffix _summary_table.txt \
      --true_tree_newick "(A,((((B,C),(D,E)),F),(G,H)));" \
      --rootfolder $rootfolder \
      --selection_method $selection_method \
      --threshold 1000`)
end
```

**Arguments accepted by `findgraphs_1rep.R`:**
- `--input_dir` (`-i`): Directory containing input sequence alignments (required)
- `--output_dir` (`-o`): Directory for output files (required)
- `--prefix` (`-p`): Input file prefix (e.g., "concatenated")
- `--num_admix` (`-k`): Number of admixture events (default: 0)
- `--stop_gen`: Number of generations to stop find_graphs (default: 100)
- `--outgroup`: Outgroup population name (default: "A")
- `--runs`: Number of independent runs (default: 100)
- `--blgsize` (`-b`): Block size used in find_graphs (default: 300)
- `--rep_id` (`-r`): Replication ID for result saving
- `--seed_file_path` (`-s`): Path to the seed file to be used
- `--output_graph_suffix`: Suffix for output graph files (default: "_graphs.rds")
- `--output_f2_suffix`: Suffix for output f2 files (default: "_f2.rds")
- `--output_summary_table_suffix`: Suffix for summary table outputs (default: "_summary_table.txt")
- `--true_tree_newick`: True species tree in Newick format (default: "(A,((((B,C),(D,E)),F),(G,H)));")
- `--rootfolder`: Root folder to save any warnings (default: "__use_cwd__")
- `--selection_method`: Method to select graphs (optional)
- `--threshold`: Threshold for graph selection (default: 1000)

### Parallel Processing

The script supports multiple levels of parallelization:
- Each replicate processed independently
- Launch with: `julia -p N findgraphs.jl --params DUP0.01-LOS0.005-RVGL-N_ind2`

### Usage Example of findgraphs.jl 

```bash
# Analyze population graphs for specific parameter set
julia -p 8 scripts/findgraphs.jl \
  --params DUP0.01-LOS0.005-RVGL-N_ind2 \
  --min_maf 0.05 \
  --test_admixture \
  --test_topology \
  --bootstrap_reps 100

# Process multiple parameter combinations
for params in DUP0.01-LOS0.005-RVGL-N_ind2 DUP0.02-LOS0.01-RVGL-N_ind2; do
    julia -p 8 scripts/findgraphs.jl --params $params
done
```

This analyzes population structure and admixture patterns in the simulated datasets, providing complementary insights to the phylogenetic network analysis from `snaq_submit.jl`.


## Stage 2.2: Network Inference (snaq_submit.jl)

### Overview

[`snaq_submit.jl`](snaq_submit.jl) uses the gene trees generated by `simulation_iqtree.jl` to infer phylogenetic networks that can capture reticulation events like hybridization and introgression. This script processes the estimated gene trees to construct species networks using the SNaQ (Species Networks applying Quartets) method implemented in PhyloNetworks.jl.

### Scripts called by snaq_submit.jl

- [`snaq_1rep.jl`](snaq_1rep.jl): Core processing script that handles individual replicates
- [`utilities.jl`](utilities.jl): Helper functions for tree processing, data management, and network analysis
- Additional utility scripts for quartet sampling and network evaluation

### Software dependencies

The pipeline uses the following Julia packages and external software:

- **[PhyloNetworks.jl](https://github.com/crsl4/PhyloNetworks.jl)**: Phylogenetic network inference using SNaQ
- **[PhyloTrees.jl](https://github.com/jangevaare/PhyloTrees.jl)**: Tree manipulation and analysis
- **[DataFrames.jl](https://github.com/JuliaData/DataFrames.jl)**: Data processing and management
- **[CSV.jl](https://github.com/JuliaData/CSV.jl)**: File I/O operations
- **[StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl)**: Statistical utilities

### Detailed Workflow

#### 1. Input Data Processing

**Purpose**: Prepare gene trees from IQ-TREE output for network inference.

**Process**:
- Reads gene trees from `iqtreefolder/` (output of `simulation_iqtree.jl`)
- Validates tree format and topology consistency
- Handles multiple individuals per species when `n_inds > 1`
- Creates species mapping files for multi-individual datasets
- Filters trees based on quality criteria (branch support, topology)

#### 2. Network Inference via snaq_1rep.jl

**Purpose**: Process individual replicates to infer phylogenetic networks with different numbers of hybridization events.

**How snaq_1rep.jl is Called**:
```julia
# For each replicate directory
for rep_dir in replicate_directories
    # Skip replicates with insufficient data (from warning logs)
    if has_warnings(rep_dir)
        continue
    end
    
    # Call snaq_1rep.jl with specific parameters
    run(`julia snaq_1rep.jl 
         --input_dir $rep_dir/iqtreefolder 
         --output_dir $rep_dir/snaq_output
         --species_mapping $mapping_file
         --max_reticulations $max_h
         --runs_per_h $n_runs
         --seed $replicate_seed`)
end
```

**snaq_1rep.jl Processing Steps**:

1. **Gene Tree Preparation**:
   - Loads gene trees from IQ-TREE output files
   - Converts tip names to species identifiers (removes individual IDs)
   - Creates quartet frequency tables for SNaQ input
   - Handles missing taxa and incomplete lineage sorting

2. **Network Inference Across Hybridization Levels**:
   - **h=0**: Infer species tree (no reticulation events)
   - **h=1**: Infer network with 1 hybridization event

3. **Multiple Independent Runs**:
   - Runs multiple SNaQ searches per hybridization level (default: 10 runs)
   - Uses different random starting points for each run
   - Implements parallel processing when multiple cores available
   - Tracks convergence and optimization statistics

4. **Network Selection and Evaluation**:
   - Compares networks using pseudo-likelihood scores
   - Selects best network per hybridization level
   - Computes network statistics (reticulation confidence, edge lengths)
   - Evaluates network stability across runs

### Input Requirements

**From simulation_iqtree.jl output**:
- Gene trees (`*.treefile`) in each replicate's `iqtreefolder/`
- Species tree inferred from wAstral (asral.tre) in `astralfolder/`
- Species/individual mapping information (`astralfolder/astral_mapping.txt`)
- Simulation parameters for result organization


### Output Structure

For each parameter combination, creates network analysis folders within existing simulation directories:

```
output/DUP{rate}-LOS{rate}-RV{var}-N_ind{n}/
├── rep001/
│   ├── iqtreefolder/           # Input: IQ-TREE gene trees
│   └── snaq_output/            # Output: SNaQ network inference
│       ├── gene_trees_input/   # Processed gene trees for SNaQ
│       ├── networks_h0/        # Species trees (h=0)
│       ├── networks_h1/        # Networks with 1 hybridization
│       ├── logs/              # SNaQ optimization logs
│       └── summary/           # Network statistics and comparisons
```
### Parallel Processing

The script supports multiple levels of parallelization:
- Each replicate processed independently
- Launch with: `julia -p N snaq_submit.jl --params DUP0.01-LOS0.005-RVGL-N_ind2`

### Seed Management

**Deterministic seeding system**:
- Uses master seed based on simulation parameters
- Each replicate gets independent seed derived from master
- SNaQ runs within replicates use sequential seeds
- Ensures reproducible network inference results

### Usage Example
```bash
# Run network inference for a specific simulation parameter set with 8 cores for snaq_submitl.jl 
julia -p 8 scripts/snaq_submit.jl \
   --params DUP0.01-LOS0.005-RVGL-N_ind2 \
   --max_reticulations 1 \
   --runs_per_h 10
```
This will infer phylogenetic networks with up to 1 hybridization event for each replicate, running 10 independent SNaQ searches per value of h, and organize results by simulation parameter set.

## snaq_submit.jl Usage Example

To run the full SNaQ network inference pipeline across all replicates for a given simulation parameter set, use the `snaq_submit.jl` script. This script manages batch processing, parallelization, and reproducibility.

```bash
# Basic usage: infer networks with up to 1 hybridization event for all replicates
julia -p 8 scripts/snaq_submit.jl \
   --params DUP0.01-LOS0.005-RVGL-N_ind2 \
   --max_reticulations 1 \
   --runs_per_h 10
```

This command will process all replicates in the specified parameter set, running 10 independent SNaQ searches for each value of h (number of hybridization events), and organize results in the appropriate output directories.


