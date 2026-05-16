# ============================================================================
# scripts/clean.jl
#
# Purpose : Remove the per-replicate output/ directory (every parameter
#           setting and replicate is wiped) plus any stray raxml-outfiles /
#           astral-outfiles folders at the repo root. Use to reset a working
#           tree before re-running the pipeline from scratch.
# Inputs  : None.
# Outputs : Removes output/ and matching folders if present.
# Usage   : julia --project=. scripts/clean.jl
# Warning : Destructive. Aggregated results in results/, simulation_summary/,
#           snaq_summary/, and findgraph_summary/ are NOT deleted, but every
#           per-replicate intermediate under output/ is.
# ============================================================================

# Remove the single output folder
rm("output", force=true, recursive=true)

# See simulation.jl, temp outputs of raxml and astral were moved later. 
# In cases if we want to remove temp outfiles from the root. 
# Check if the folder exists and then remove them. 
# If not exist, this step is uncessary. 
names_to_be_removed = ["raxml-outfiles", "astral-outfiles"]
for folder in readdir()
    if any(occursin(name, folder) for name in names_to_be_removed)
        rm(folder; force = true, recursive = true)
        println("Removed temp files from $folder in the root folder")
    end
end