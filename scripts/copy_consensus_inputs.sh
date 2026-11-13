#!/bin/bash

# Function of this script: 
# copy the few network files needed for the consensus trees / networks
# (snaq_postprocess.jl and findgraphs_postprocess.jl) from the servers to here,
# keeping the folder structure. Run this first, then the postprocess scripts.
# The parameter folders (DUP*-LOS*-RV*-N_ind*-SF*-genelen*) 
# are found in the input folder. 

# Thus, it was mainly used for transferring files, 
# and *not* used in the main pipeline. 

# Usage: bash copy_consensus_inputs.sh --mode <snaq|findgraphs> [options], see usage() below
# snaq copies, for each rep: snaqfolder/H0_output/H0.out and snaqfolder/H1_output/H1.out
# findgraphs copies, for each rep: findgraph/rep<id>_admix{0,1}_unique_graphs.rds,
#   findgraph/rep<id>_admix{0,1}_summary_table.txt, findgraph/rep<id>_f2.rds
#   and once per parameter set: findgraphs_summary_results.csv

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") --mode <snaq|findgraphs> [options]

Required:
  --mode snaq | findgraphs

Optional (source / destination):
  --input_folder folders that contain the DUP*-LOS*-RV*... parameter folders,
      several folders separated by : (default: <dest>/output). --sources is an old name for it
  --dest destination root (default: the repo root, one level up from this script)

Optional (parameter filters, if omitted all parameter sets found are processed):
  --dup_rate duplication rate (e.g. 0.0003)
  --loss_rate gene loss rate (e.g. 0.0003)
  --ratevar rate variation as in the folder name, after RV: N, L, G or G*L
  --n_inds individuals per species (1 or 2)
  --SF scaling factor (e.g. 1.0)
  --gene_len gene length (e.g. 1000)

Optional (replicate range):
  --rep_start first replicate to copy (default: all)
  --rep_end last replicate to copy (default: all)

Other:
  --dry_run print what would be copied, copy nothing
  -h | --help show this message

Examples:
  # Copy everything found under the default output/ folder:
  bash copy_consensus_inputs.sh --mode findgraphs

  # Copy only the parameter sets without rate variation (RVN) from a remote source path:
  bash copy_consensus_inputs.sh --mode snaq --ratevar N \\
      --input_folder /path/to/remote/sim/output

  # Multiple source roots (colon-separated), with a dup_rate filter:
  bash copy_consensus_inputs.sh --mode findgraphs --dup_rate 0.0003 \\
      --input_folder "/path/to/server1/output:/path/to/server2/output"
EOF
    exit 1
}

# defaults
MODE=""
INPUT_FOLDER="" # set below after DEST, if still empty
DUP_RATE=""
LOSS_RATE=""
RATEVAR=""
N_INDS=""
SF=""
GENE_LEN=""
REP_START=""
REP_END=""
# default DEST = repo root, one level up from scripts/. --dest to change it
DEST="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DRY_RUN=false

# arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --input_folder) INPUT_FOLDER="$2"; shift 2 ;;
        --sources) INPUT_FOLDER="$2"; shift 2 ;; # old name of --input_folder
        --dup_rate) DUP_RATE="$2"; shift 2 ;;
        --loss_rate) LOSS_RATE="$2"; shift 2 ;;
        --ratevar) RATEVAR="$2"; shift 2 ;;
        --n_inds) N_INDS="$2"; shift 2 ;;
        --SF) SF="$2"; shift 2 ;;
        --gene_len) GENE_LEN="$2"; shift 2 ;;
        --rep_start) REP_START="$2"; shift 2 ;;
        --rep_end) REP_END="$2"; shift 2 ;;
        --dest) DEST="$2"; shift 2 ;;
        --dry_run) DRY_RUN=true; shift 1 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# check arguments
if [[ -z "$MODE" ]]; then
    echo "Error: --mode is required"
    usage
fi

if [[ "$MODE" != "snaq" && "$MODE" != "findgraphs" ]]; then
    echo "Error: --mode must be 'snaq' or 'findgraphs' (got: '$MODE')"
    usage
fi

# default input folder: output/ in DEST
[[ -z "$INPUT_FOLDER" ]] && INPUT_FOLDER="$DEST/output"

# copy one file (src_file) to $DEST/output/<path relative to input_root>
# cp -p keeps the timestamps, -u skips files where the destination is newer
do_copy() {
    local src_file="$1"
    local input_root="$2"
    local rel_path="${src_file#${input_root}/}"
    local dest_file="$DEST/output/$rel_path"
    local dest_dir
    dest_dir=$(dirname "$dest_file")

    # parameter folder already at the destination: 
    # only copy newer files (-u), otherwise copy everything
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

# several input roots can be given, separated by :
IFS=':' read -ra INPUT_ROOTS <<< "$INPUT_FOLDER"

# does this parameter folder name pass all the filters given?
param_matches_filters() {
    local pname="$1"
    [[ -n "$DUP_RATE" && "$pname" != *"DUP${DUP_RATE}-"* ]] && return 1
    [[ -n "$LOSS_RATE" && "$pname" != *"LOS${LOSS_RATE}-"* ]] && return 1
    [[ -n "$RATEVAR" && "$pname" != *"-RV${RATEVAR}-"* ]] && return 1
    [[ -n "$N_INDS" && "$pname" != *"N_ind${N_INDS}-"* ]] && return 1
    [[ -n "$SF" && "$pname" != *"-SF${SF}-"* ]] && return 1
    [[ -n "$GENE_LEN" && "$pname" != *"genelen${GENE_LEN}"* ]] && return 1
    return 0
}

# all parameter set names found in the input roots
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

echo "mode: $MODE"
echo "input folder(s): $INPUT_FOLDER"
echo "destination: $DEST"
echo "dry run: $DRY_RUN"
echo "parameter sets: ${#PARAM_SEEN[@]}"
[[ -n "$DUP_RATE" ]] && echo "filter dup_rate: $DUP_RATE"
[[ -n "$LOSS_RATE" ]] && echo "filter loss_rate: $LOSS_RATE"
[[ -n "$RATEVAR" ]] && echo "filter ratevar: $RATEVAR"
[[ -n "$N_INDS" ]] && echo "filter n_inds: $N_INDS"
[[ -n "$SF" ]] && echo "filter SF: $SF"
[[ -n "$GENE_LEN" ]] && echo "filter gene_len: $GENE_LEN"
[[ -n "$REP_START" || -n "$REP_END" ]] && \
    echo "rep range: ${REP_START:-*} to ${REP_END:-*}"

n_copied=0
n_missing=0
n_reps_not_found=0

# go through the parameter sets
for pname in $(printf '%s\n' "${!PARAM_SEEN[@]}" | sort); do
    echo ""
    echo "--- $pname ---"

    # rep folders in all input roots. If a rep is in several roots,
    # the first one wins
    declare -A REP_ROOT=()
    for input_root in "${INPUT_ROOTS[@]}"; do
        param_path="$input_root/$pname"
        [[ -d "$param_path" ]] || continue
        for rep_dir in "$param_path"/rep*; do
            [[ -d "$rep_dir" ]] || continue
            rep_id=$(basename "$rep_dir")
            [[ -z "${REP_ROOT[$rep_id]:-}" ]] && \
                REP_ROOT["$rep_id"]="$input_root"
        done
    done

    if [[ ${#REP_ROOT[@]} -eq 0 ]]; then
        echo "  WARNING: no replicate directories found, skipping"
        continue
    fi

    # copy the files of each rep
    for rep_id in $(printf '%s\n' "${!REP_ROOT[@]}" | sort); do
        # rep range filter: rep number without "rep" and without leading zeros
        rep_num="${rep_id#rep}"
        rep_num="${rep_num#"${rep_num%%[!0]*}"}"
        rep_num=$((10#${rep_num}))
        [[ -n "$REP_START" && "$rep_num" -lt "$REP_START" ]] && continue
        [[ -n "$REP_END" && "$rep_num" -gt "$REP_END" ]] && continue

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

    # findgraphs only: 
    # findgraphs_summary_results.csv at the top of the parameter folder,
    # written by findgraphs_1rep.R and read by generate_summary_csv 
    # in findgraphs_postprocess.jl
    if [[ "$MODE" == "findgraphs" ]]; then
        toplevel_files=( "findgraphs_summary_results.csv" )
        for fname in "${toplevel_files[@]}"; do
            for input_root in "${INPUT_ROOTS[@]}"; do
                candidate="$input_root/$pname/$fname"
                if [[ -f "$candidate" ]]; then
                    echo "top-level: $fname  [source: $input_root]"
                    do_copy "$candidate" "$input_root"
                    (( n_copied++ )) || true
                    break   # first root that has it
                fi
            done
        done
    fi

    unset REP_ROOT
    declare -A REP_ROOT=()
done

# summary
echo "done. files copied or already up to date: $n_copied"
echo "files missing on the source: $n_missing"
echo "replicate directories not found: $n_reps_not_found"
