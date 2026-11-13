# Remove output/ (all parameter settings and replicates) and the leftover
# raxml-outfiles / astral-outfiles folders at the root, to start over.
# Usage: julia scripts/clean.jl
# Careful: everything under output/ is deleted. The summaries in results/,
# simulation_summary/, snaq_summary/ and findgraph_summary/ are kept

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