# Executables

This directory holds to the third-party binaries required to run the simulation pipeline. The actual binaries are installed separately (see [`software_installation.sh`](software_installation.sh)) and linked here so that all pipeline scripts can reference them through a single, stable relative path (`executables/<name>`).

> **Note:** `Project.toml` and `Manifest.toml` are at the **repository root**

---

## Required Executables

| Symlink | Software | Version | Purpose |
|---|---|---|---|
| `simphy` | [SimPhy](https://github.com/adamallo/SimPhy) | v1.0.2 | Simulate gene trees along a species phylogeny (with ILS, duplication, loss) |
| `seq-gen` | [Seq-Gen](https://github.com/rambaut/Seq-Gen) | v1.3.5 | Simulate DNA sequences along gene trees |
| `iqtree2` | [IQ-TREE 2](https://github.com/iqtree/iqtree2) | v2.4.0 | Maximum-likelihood gene tree estimation from simulated sequences |
| `astral-IV` | [ASTER](https://github.com/chaoszhang/ASTER) | v1.24.4.8 | Species tree estimation from gene trees (ASTRAL algorithm) |
| `snp-sites` | [snp-sites](https://github.com/sanger-pathogens/snp-sites) | v2.5.1 | Extract SNPs from whole-genome alignments produced by Seq-Gen |

*Note*: We used wastral and astral-pro for testing only. 

All entries in this directory are **symbolic links** to binaries installed in a private software directory (e.g., `~/private/software/`). The binaries themselves are **not tracked by Git** (see [`.gitignore`](../.gitignore)).

---

## Installation

Detailed, step-by-step installation instructions for each binary are provided in [`software_installation.sh`](software_installation.sh). The general pattern for each tool is:

1. Download and build/unpack the release archive into a private software directory.
2. Create a symlink inside this `executables/` directory pointing to the installed binary.
3. Verify the installation by running the binary with a help or version flag.

```bash
# General pattern (replace <tool> and <path> accordingly)
cd ~/private/software/
# ... download and build <tool> ...
cd /path/to/simulation-reptiles
ln -s ~/private/software/<tool-install-path>/<binary> executables/<symlink-name>
./executables/<symlink-name> --help # Just to verify if the link works 
```

For the full commands including exact download URLs, version pins, and build flags, see [`software_installation.sh`](software_installation.sh).

---

## Platform

All binaries were compiled and tested on **Linux x86-64** (UW–Madison statistics department servers, e.g., `franklin00`, `franklin01`, `franklin02`, `franklin03`).
They are not expected to run on macOS or Windows without recompilation.

---

## Additional Runtime Dependencies

Beyond the binaries above, the pipeline also requires:

| Dependency | Version | Install method |
|---|---|---|
| **Julia** | ≥ 1.11 | `curl -fsSL https://install.julialang.org \| sh` |
| Julia packages | pinned in `../Project.toml` / `../Manifest.toml` | `julia --project=.. -e 'using Pkg; Pkg.instantiate()'` |
| **Python** | 3.12 | Built from source (see `software_installation.sh`) |
| `biopython` | latest | `pip install biopython` |
| **R** | system R | — |
| `admixtools` | 2.0.8 | `devtools::install_github("uqrmaie1/admixtools")` |
| `dplyr` | 1.1.4 | `install.packages("dplyr")` |
| `optparse` | 1.7.5 | `install.packages("optparse")` |
| `Rcpp` | 1.0.14 | `install.packages("Rcpp")` |

See [`software_installation.sh`](software_installation.sh) for exact installation commands.
