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
        println("Removed temporary files from $folder stored in the root folder")
    end
end