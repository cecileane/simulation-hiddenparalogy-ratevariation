#!/bin/bash
# copy_consensus_inputs.sh
#
# Copy only the minimal set of network output files needed for consensus tree /
# network computation from a source folder to a destination, preserving the
# original directory structure.
#
# Parameter combinations (DUP*/LOS*/RV*/N_ind*/SF*/genelen*) are discovered
# automatically from the input folder — no need to spell out dup_rate,
# loss_rate, ratevar, etc.  Optional filter flags narrow which combinations
# are processed.
#
# Designed to work alongside snaq_postprocess.jl and findgraphs_postprocess.jl.
# Run this script first to stage the necessary files, then run the Julia
# postprocessing scripts on the destination as normal.
#
# Usage:
#   bash copy_consensus_inputs.sh --mode <snaq|findgraphs> [options]
#
# Files copied per mode
# ----------------------
# snaq:
#   output/<paramname_root>/rep<id>/snaqfolder/H0_output/H0.out
#   output/<paramname_root>/rep<id>/snaqfolder/H1_output/H1.out
#
# findgraphs:
#   output/<paramname_root>/rep<id>/findgraph/rep<id>_admix0_unique_graphs.rds
#   output/<paramname_root>/rep<id>/findgraph/rep<id>_admix0_summary_table.txt
#   output/<paramname_root>/rep<id>/findgraph/rep<id>_admix1_unique_graphs.rds
#   output/<paramname_root>/rep<id>/findgraph/rep<id>_admix1_summary_table.txt
#   output/<paramname_root>/rep<id>/findgraph/rep<id>_f2.rds
#   output/<paramname_root>/findgraphs_summary_results.csv  (top-level, once)

set -euo pipefail

# ---------------------------------------------------------------------------
# Help / usage
# ---------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: $(basename "$0") --mode <snaq|findgraphs> [options]

Required:
  --mode          snaq | findgraphs

Optional (source / destination):
  --input_folder  Colon-separated list of folders that directly contain the
                  DUP*/LOS*/RV*... parameter subdirectories.
                  Default: <dest>/output
                  Alias: --sources (kept for backward compatibility)
  --dest          Destination root path
                  (default: /u/b/i/bingl/private/simulation-reptiles)

Optional (parameter filters — if omitted, all discovered values are processed):
  --dup_rate      Filter by duplication rate        (e.g. 0.0003)
  --loss_rate     Filter by gene loss rate          (e.g. 0.0003)
  --ratevar       Filter by rate variation label    (e.g. RVN, RVL, RVG)
  --n_inds        Filter by individuals per species (e.g. 1)
  --SF            Filter by branch-length scale     (e.g. 1.0)
  --gene_len      Filter by gene length             (e.g. 1000)

Optional (replicate range):
  --rep_start     First replicate to copy  (default: process all found)
  --rep_end       Last  replicate to copy  (default: process all found)

Other:
  --dry_run       Print what would be copied without actually copying
  -h | --help     Show this message and exit

Examples:
  # Copy everything found under the default output/ folder:
  bash copy_consensus_inputs.sh --mode findgraphs

  # Copy only RVN parameter sets from a remote path:
  bash copy_consensus_inputs.sh --mode snaq --ratevar RVN \\
      --input_folder /nobackup/bingli/simulation-reptiles/output

  # Multiple source roots (colon-separated), with a dup_rate filter:
  bash copy_consensus_inputs.sh --mode findgraphs --dup_rate 0.0003 \\
      --input_folder "/nobackup/bingli/sim/output:/nobackup2/bingli/sim/output"
EOF
    exit 1
}

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
MODE=""
INPUT_FOLDER=""   # resolved after DEST is set, if still empty
DUP_RATE=""
LOSS_RATE=""
RATEVAR=""
N_INDS=""
SF=""
GENE_LEN=""
REP_START=""
REP_END=""
DEST="/u/b/i/bingl/private/simulation-reptiles"
DRY_RUN=false

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)           MODE="$2";         shift 2 ;;
        --input_folder)   INPUT_FOLDER="$2"; shift 2 ;;
        --sources)        INPUT_FOLDER="$2"; shift 2 ;;  # backward-compat alias
        --dup_rate)       DUP_RATE="$2";     shift 2 ;;
        --loss_rate)      LOSS_RATE="$2";    shift 2 ;;
        --ratevar)        RATEVAR="$2";      shift 2 ;;
        --n_inds)         N_INDS="$2";       shift 2 ;;
        --SF)             SF="$2";           shift 2 ;;
        --gene_len)       GENE_LEN="$2";     shift 2 ;;
        --rep_start)      REP_START="$2";    shift 2 ;;
        --rep_end)        REP_END="$2";      shift 2 ;;
        --dest)           DEST="$2";         shift 2 ;;
        --dry_run)        DRY_RUN=true;      shift 1 ;;
        -h|--help)        usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate required arguments
# ---------------------------------------------------------------------------
if [[ -z "$MODE" ]]; then
    echo "Error: --mode is required"
    usage
fi

if [[ "$MODE" != "snaq" && "$MODE" != "findgraphs" ]]; then
    echo "Error: --mode must be 'snaq' or 'findgraphs' (got: '$MODE')"
    usage
fi

# Default input folder: output/ subdirectory inside DEST
[[ -z "$INPUT_FOLDER" ]] && INPUT_FOLDER="$DEST/output"

# ---------------------------------------------------------------------------
# Helper: copy one file, preserving relative structure under dest root.
#   src_file    : absolute path on the source side
#   input_root  : the input folder root to strip when building the dest path
#                 (files are placed under $DEST/<rel_from_input_root>)
# Uses `cp -pu` (preserve timestamps, skip if dest is newer).
# ---------------------------------------------------------------------------
do_copy() {
    local src_file="$1"
    local input_root="$2"
    local rel_path="${src_file#${input_root}/}"
    local dest_file="$DEST/output/$rel_path"
    local dest_dir
    dest_dir=$(dirname "$dest_file")

    # If the parameter set folder already exists at the destination, use -u
    # (skip files where the destination is already newer).  Otherwise copy
    # unconditionally so the whole folder is transferred fresh.
    local param_dir_name="${rel_path%%/*}"
    local cp_flags="-p"
    [[ -d "$DEST/output/$param_dir_name" ]] && cp_flags="-pu"

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY RUN] output/$rel_path  (cp flags: $cp_flags)"
        echo "            -> $dest_file"
    else
        mkdir -p "$dest_dir"
        if cp $cp_flags "$src_file" "$dest_file"; then
            echo "  OK  output/$rel_path"
        else
            echo "  FAIL output/$rel_path"
            return 1
        fi
    fi
}

# ---------------------------------------------------------------------------
# Split INPUT_FOLDER on ':' into an array of source roots
# ---------------------------------------------------------------------------
IFS=':' read -ra INPUT_ROOTS <<< "$INPUT_FOLDER"

# ---------------------------------------------------------------------------
# Helper: check whether a param directory name matches all active filters
# ---------------------------------------------------------------------------
param_matches_filters() {
    local pname="$1"
    [[ -n "$DUP_RATE"  && "$pname" != *"DUP${DUP_RATE}-"*   ]] && return 1
    [[ -n "$LOSS_RATE" && "$pname" != *"LOS${LOSS_RATE}-"*  ]] && return 1
    [[ -n "$RATEVAR"   && "$pname" != *"-RV${RATEVAR}-"*    ]] && return 1
    [[ -n "$N_INDS"    && "$pname" != *"N_ind${N_INDS}-"*   ]] && return 1
    [[ -n "$SF"        && "$pname" != *"-SF${SF}-"*          ]] && return 1
    [[ -n "$GENE_LEN"  && "$pname" != *"genelen${GENE_LEN}"* ]] && return 1
    return 0
}

# ---------------------------------------------------------------------------
# Discover all unique parameter-set names across all input roots
# ---------------------------------------------------------------------------
declare -A PARAM_SEEN=()
for input_root in "${INPUT_ROOTS[@]}"; do
    if [[ ! -d "$input_root" ]]; then
        echo "WARNING: input folder not found, skipping: $input_root"
        continue
    fi
    for param_dir in "$input_root"/DUP*-LOS*-RV*-N_ind*-SF*-genelen*; do
        [[ -d "$param_dir" ]] || continue
        pname=$(basename "$param_dir")
        param_matches_filters "$pname" && PARAM_SEEN["$pname"]=1
    done
done

if [[ ${#PARAM_SEEN[@]} -eq 0 ]]; then
    echo "No parameter directories found matching the given filters in:"
    printf '  %s\n' "${INPUT_ROOTS[@]}"
    exit 1
fi

echo "============================================"
echo "  Mode              : $MODE"
echo "  Input folder(s)   : $INPUT_FOLDER"
echo "  Destination root  : $DEST"
echo "  Dry run           : $DRY_RUN"
echo "  Parameter sets    : ${#PARAM_SEEN[@]}"
[[ -n "$DUP_RATE"  ]] && echo "  Filter dup_rate   : $DUP_RATE"
[[ -n "$LOSS_RATE" ]] && echo "  Filter loss_rate  : $LOSS_RATE"
[[ -n "$RATEVAR"   ]] && echo "  Filter ratevar    : $RATEVAR"
[[ -n "$N_INDS"    ]] && echo "  Filter n_inds     : $N_INDS"
[[ -n "$SF"        ]] && echo "  Filter SF         : $SF"
[[ -n "$GENE_LEN"  ]] && echo "  Filter gene_len   : $GENE_LEN"
[[ -n "$REP_START" || -n "$REP_END" ]] && \
    echo "  Rep range         : ${REP_START:-*} to ${REP_END:-*}"
echo "============================================"

n_copied=0
n_missing=0
n_reps_not_found=0

# ---------------------------------------------------------------------------
# Process each discovered parameter set
# ---------------------------------------------------------------------------
for pname in $(printf '%s\n' "${!PARAM_SEEN[@]}" | sort); do
    echo ""
    echo "--- $pname ---"

    # Collect all rep directories across input roots (first root wins per rep)
    declare -A REP_ROOT=()
    for input_root in "${INPUT_ROOTS[@]}"; do
        param_path="$input_root/$pname"
        [[ -d "$param_path" ]] || continue
        for rep_dir in "$param_path"/rep*; do
            [[ -d "$rep_dir" ]] || continue
            rep_id=$(basename "$rep_dir")
            # first source root that has this rep wins
            [[ -z "${REP_ROOT[$rep_id]:-}" ]] && \
                REP_ROOT["$rep_id"]="$input_root"
        done
    done

    if [[ ${#REP_ROOT[@]} -eq 0 ]]; then
        echo "  WARNING: no replicate directories found — skipping"
        continue
    fi

    # -----------------------------------------------------------------
    # Per-replicate file copy
    # -----------------------------------------------------------------
    for rep_id in $(printf '%s\n' "${!REP_ROOT[@]}" | sort); do
        # Apply optional rep range filters
        rep_num="${rep_id#rep}"   # strip leading "rep" to get numeric part
        rep_num="${rep_num#"${rep_num%%[!0]*}"}"  # strip leading zeros
        rep_num=$((10#${rep_num}))                # force base-10
        [[ -n "$REP_START" && "$rep_num" -lt "$REP_START" ]] && continue
        [[ -n "$REP_END"   && "$rep_num" -gt "$REP_END"   ]] && continue

        input_root="${REP_ROOT[$rep_id]}"
        src_rep_dir="$input_root/$pname/$rep_id"

        echo "$rep_id  [source: $input_root]"

        declare -a files_to_copy=()

        if [[ "$MODE" == "snaq" ]]; then
            files_to_copy=(
                "$src_rep_dir/snaqfolder/H0_output/H0.out"
                "$src_rep_dir/snaqfolder/H1_output/H1.out"
            )
        elif [[ "$MODE" == "findgraphs" ]]; then
            files_to_copy=(
                "$src_rep_dir/findgraph/${rep_id}_admix0_unique_graphs.rds"
                "$src_rep_dir/findgraph/${rep_id}_admix0_summary_table.txt"
                "$src_rep_dir/findgraph/${rep_id}_admix1_unique_graphs.rds"
                "$src_rep_dir/findgraph/${rep_id}_admix1_summary_table.txt"
                "$src_rep_dir/findgraph/${rep_id}_f2.rds"
            )
        fi

        for f in "${files_to_copy[@]}"; do
            if [[ -f "$f" ]]; then
                do_copy "$f" "$input_root"
                (( n_copied++ )) || true
            else
                echo "  MISSING  $f"
                (( n_missing++ )) || true
            fi
        done
    done

    # -------------------------------------------------------------------------
    # Top-level files (findgraphs only)
    # findgraphs_summary_results.csv is generated by the R findgraphs pipeline
    # and is read by generate_summary_csv() in findgraphs_postprocess.jl
    # -------------------------------------------------------------------------
    if [[ "$MODE" == "findgraphs" ]]; then
        toplevel_files=( "findgraphs_summary_results.csv" )
        for fname in "${toplevel_files[@]}"; do
            for input_root in "${INPUT_ROOTS[@]}"; do
                candidate="$input_root/$pname/$fname"
                if [[ -f "$candidate" ]]; then
                    echo "top-level: $fname  [source: $input_root]"
                    do_copy "$candidate" "$input_root"
                    (( n_copied++ )) || true
                    break   # only copy from the first source that has it
                fi
            done
        done
    fi

    unset REP_ROOT
    declare -A REP_ROOT=()
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "============================================"
echo "  Done."
echo "  Files copied / already up-to-date : $n_copied"
echo "  Expected files missing on source  : $n_missing"
echo "  Replicate directories not found   : $n_reps_not_found"
echo "============================================"
