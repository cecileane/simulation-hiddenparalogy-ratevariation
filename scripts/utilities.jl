#= Script to store utility functions 
The below script has four parts 

Part 0 --> Logging helper functions for distributed computation
    F1: "worker_println" 
        -> Print to stdout and log file if logging is enabled
    F2: "open_worker_log" 
        -> Open worker-specific log file
    F3: "close_worker_log" 
        -> Close worker-specific log file

Part 1 --> A group of functions to modify gene tree output from SimPhy
    F1: "replace_taxa_name" 
        -> Count and modify one tree tip given Char and newick 
    F2: "count_missing" 
        -> Count # missing taxa in the given newick and modify 
        all tree tips (conditions seen below)
    F3: "modify_newicks" 
        -> write modified newicks (condtions seen below)
    F4: "modify_newick_for_n_genes" 
        -> run "modify_newicks" for n times and output the files with given prefix 

Part 2 --> A group of helper functions to re-run simphy simulation and modify newicks
 PS:       seen details in scripts/speciestree_iqtree.jl 
    F1: “seed_generator” 
        -> generate N x M arrays with random seeds given a master seed 
    F2: "calculate_batch" 
        -> Estimate batch size for simulating gene trees in the next iterations 

Part 3 --> Other utilities functions 
    F1: "pad_number" 
        -> convert a given Int to a string to meet the output from Simphy
    F2: "check_existing_dir" 
        -> check if a path exists and then decide to remove the path or not
    F3: "replace_tips_with_letters" 
        -> Replace text tips with letters to both strings and Network/tree 
        seen scripts/speciestree.jl 
    F4: "replace_pop_inInd_files" 
        -> Replace "POP" from the output ind file to make each sample is 
        mapped to a unique population ID
    F5: "generate_software_seeds" 
        -> Generate a dictionary of stable random seeds for each software 
        name using a given master seed.
    F6: "generate_master_seed" 
        -> Generate a master seed based on given parameter values.
    F7: "get_dict_for_seed_setting" 
        -> Get a dictionary for setting up the master seed based on the 
        string "paramname_root"
    F8: "median_branch_length" 
        -> Compute the median branch length of a phylogenetic tree or network.
    F9: "transform_newick_no_lineage_tree" 
        -> Transform a Newick tree string by scaling branch lengths based on 
        a given pattern and effective population size (Ne). 
    F10: "transform_newick_with_lineage_tree" 
        -> Transform a Newick tree string by scaling branch lengths based on 
        lineage tree data.
    F11: "add_outgroup_to_newick" 
        -> Add an outgroup branch to a Newick tree string. 
    F12: "prune_tip_from_newick" 
        -> Remove a specific tip from a Newick tree string. 
    F13: "totallength" 
        -> Calculate the total branch length of a phylogenetic tree or network.  
    F14: "hasrep" 
        -> Check if a Newick string has repeated taxa. 
        This helper function is retrieved from SNaQ.jl documentation 
            https://juliaphylo.github.io/SNaQ.jl/dev/man/multiplealleles/ 
    F15: "set_up_paramname_root"
        -> Set up the string "paramname_root" based on given parameter values.
        # Super important for almost all scripts in this project 
    F16: "scale_branch_lengths"
        -> Scale the branch lengths of a Newick tree string by a given factor.
        # Used in scripts/simulation.jl
=# 
using PhyloNetworks # Used in "replace_tips_with_letters" to get trees and networks 
using StatsBase # Use the countmap function to find duplicated items
using StableRNGs # Used to generate stanble random seeds that won't change between julia versions
using Primes # use "nextprimes()" function to get the next prime bigger than a number 

using PhyloNetworks 
using StatsBase 
using StableRNGs 
using Primes 
using Statistics
using Printf 
using InlineStrings # use String3 
using IterTools # produce cartesian product of multiple vectors 
using Glob 
using Distributions 

#-----------------------------------------------#       
#               Part 0  
#    Logging helper functions
#-----------------------------------------------#
"""
    worker_println(args...)

Print to stdout and flush immediately. If logging is enabled (via global 
`log_output` flag and `worker_logfile`), also write to the worker-specific 
log file.

This function is designed for use in distributed worker processes to ensure
output is visible and logged.

# Arguments
- `args...`: Arguments to pass to println

# Notes
- Requires global variables: `log_output` (Bool) and `worker_logfile` (IOStream or nothing)
- Silently fails if log write fails to prevent interrupting the main computation

# Examples
```julia
worker_println("Rep1 Iter5: Processing...")
worker_println("Status: ", 10, " trees generated")
```
"""
function worker_println(args...)
    println(args...)
    flush(stdout)
    if @isdefined(log_output) && log_output && @isdefined(worker_logfile) && worker_logfile !== nothing
        try
            println(worker_logfile, args...)
            flush(worker_logfile)
        catch e
            # Silent fail if log write fails
        end
    end
end

"""
    open_worker_log()

Open a worker-specific log file for the current worker process.
Creates a log file named `screen_<paramname>_worker<id>.log` in the output folder.

This function should be called at the start of distributed computation when 
logging is enabled.

# Notes
- Requires global variables: `log_output`, `outfolder`, `paramname_root`
- Uses `myid()` to create unique log files per worker
- Warns but does not fail if log file cannot be opened

# Examples
```julia
@everywhere open_worker_log()
```
"""
function open_worker_log()
    if @isdefined(log_output) && log_output
        try
            log_filename = joinpath(outfolder, "screen_$(paramname_root)_worker$(myid()).log")
            global worker_logfile = open(log_filename, "w")
            println(worker_logfile, "=== Worker $(myid()) Log Started ===")
            flush(worker_logfile)
        catch e
            @warn "Could not open worker log file: $e"
        end
    end
end

"""
    close_worker_log()

Close the worker-specific log file for the current worker process.

This function should be called at the end of distributed computation when 
logging is enabled.

# Notes
- Requires global variable: `worker_logfile`
- Silently fails if close operation fails

# Examples
```julia
@everywhere close_worker_log()
```
"""
function close_worker_log()
    if @isdefined(worker_logfile) && worker_logfile !== nothing
        try
            println(worker_logfile, "=== Worker $(myid()) Log Ended ===")
            close(worker_logfile)
            global worker_logfile = nothing
        catch e
            # Silent fail
        end
    end
end

#-----------------------------------------------#       
#               Part 1  
#    Modify newick strings generated from Simphy
#    Goal: Mimic hidden paralogy    
#-----------------------------------------------#
"""
    replace_taxa_name(taxon::Char, newick::String)

Input:
- `char`: target character to search for in tip labels.
- `tree_str`: newick tree string containing tip labels in the format `char_x_y`.

Output: Tuple (count, modified_tree):
- `count`: number of occurrences of `char` in the tree.
- `modified_tree`: the newick string, modified if applicable.
   * If tip appears in the string only once, then the tree is modified.
   * If tip appears in the string 0 or more than once, then tree is not modified.

Examples: 

Case 1: the character appears in the string for only once
```julia
julia> replace_taxa_name('A', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
(1, "(A_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
```

Case 2: the character appears in the string for 0 time
```julia
julia> replace_taxa_name('F', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
(0, "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
```

Case 3: the character appears in the string >= 2 times:
```julia
julia> replace_taxa_name('A', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_1_0:0.02,A_2_1:0.02):0.01)));")
(2, "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_1_0:0.02,A_2_1:0.02):0.01)));")
```

Case 4: Similar to above, the character appears in the string >= 2 times,
but the there are different individuals from the same species:
```julia
julia> replace_taxa_name('A', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_0_0:0.02,A_0_2:0.02):0.01)));")
(1, "(A_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_0:0.02,A_2:0.02):0.01)));")
```
becase all three As are the same gene copies from different accessions. 
"""
function replace_taxa_name(taxon::Char, newick::String)
    tre = readnewick(newick)
    tips = tiplabels(tre)

    # replace "A_(\d+)_(\d+)" with "A_\2"
    # From simphy tutorial, tips are labelled as taxa_locusID_accessionID, 
    # here, remove the locusID while keeping the accessionID. 
    regex = Regex("($(taxon))_\\d+_(\\d+)")
    matching_tips = filter(tip -> occursin(regex, tip), tips)

    # Needs to handle if there are multiple accesions in one species 
    # eg 1: A_1 and A_0 are different accessions for one species, so keep both 
    # eg 2: A_0 and A_0 are the same accession for one species 
    #       that the genes got duplicated
    tip_list = [replace(tip, regex => s"\1_\2") for tip in matching_tips]
    counts = values(countmap(tip_list)) 
    num = isempty(counts) ? 0 : maximum(counts) 
    # The above code makes sure that it only takes the maximum. 
    # Thus, if one individual has more than one copies, 
    # the tree will be discarded later. If the char is not in the string, num = 0 
   

    if num == 0 || num >= 2 # If # of taxons >= 2 and == 0, no change to the tips: 
        return (num, newick)
    else # If # of taxons == 1, remove locus ID from the tip 
        modified_newick = replace(newick, regex => s"\1_\2") 
        return(num, modified_newick)
    end
end

"""
    count_missing(newick::String)

Check a Newick tree string for missing or duplicate taxa.

Each taxon label must follow the format `char_x_y`, where `char` is a species code 
from "A" to "H".  
The function identifies missing species and simplifies taxon labels to `char_index`.

Returns a tuple `(missing_count, modified_tree)` if all gene copies are unique.  
Returns `nothing` if any gene copy is duplicated.

# Arguments
- `newick::String`: A Newick-formatted tree with taxa in the format `char_x_y`.

# Returns
- `Tuple{Int, String}`: Num of missing taxa and the simplified Newick, if valid.
- `nothing`: If any taxon appears more than once with the same gene copy.

# Examples
```julia-repl
julia> count_missing("(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
(3, "(A_1:0.01,(((B_1:0.004,C_1:0.003):0.007,(D_0:0.02,E_1:0.02):0.01)));")

julia> count_missing("(F_0_1:0.8,(A_0_1:0.01,((A_1_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
nothing
```
"""
function count_missing(newick::String)
    taxa = "ABCDEFGH"  # Species tree only have those taxa (8 species)
    num_missing = 0  
    nwk = newick

    for taxon in taxa
        (ncopy, nwk) = replace_taxa_name(taxon, nwk)
        if ncopy > 1 # If any taxa appears for more than once, return nothing  
            return nothing  
        elseif ncopy == 0 # If missing any taxa, count the number of missing 
            num_missing += 1
        end
    end
    return (num_missing, nwk) # Return (number of missing, modified newick string) 
end

"""
    modify_newicks(input::String, output::String) -> Tuple{Int, Int}

Process a file of Newick trees to exclude trees with repeated taxa or excessive missing taxa,
and write the simplified valid trees to an output file.

Each tree is checked for two conditions:
- If it contains repeated gene copies (i.e. paralogy is not hidden), it is excluded.
- If more than four of the eight expected taxa (A–H) are missing, it is excluded.

Trees that pass both checks are modified by simplifying tip labels from the format
`A_locusID_accessionID` to `A_accessionID`, and are written to the output file.

# Arguments
- `input::String`: Path to the input file containing Newick trees (one per line).
- `output::String`: Path to the output file for valid and simplified Newick trees.

# Returns
- `Tuple{Int, Int}`: Number of excluded trees due to:
    1. Repeated taxa (non-hidden paralogy),
    2. Too many missing taxa.

# Notes
- If no hidden paralogy is simulated (`dup_rate = 0`, `loss_rate = 0`),
  no trees are excluded for duplicates.
- Valid trees will always be simplified even when no trees are excluded.

# Examples
```julia-repl
julia> modify_newicks("simphy_output.txt", "cleaned_output.txt")
(5, 12)  # 5 trees excluded due to repeated taxa, 12 due to too many missing
```
"""
function modify_newicks(input::String, output::String, 
                        intermediate_file_mapsl::String, 
                        intermediate_file_maplg::String,
                        max_taxa_missing::Int,
                        debug_mode::Bool=false)
    # Check if input file exists before attempting to read it
    if !isfile(input)
        @warn "Input file $input does not exist, skipping"
        return (0, 0)
    end
    
    # Read the entire file content atomically before any processing
    # This prevents file system race conditions on network file systems (like AFS)
    # where file operations can have delays between open/read/close
    # Note: SimPhy gene tree files contain a single Newick string per file
    nwk = try
        # Read entire file, remove whitespace, convert to String
        content = read(input, String)
        String(strip(content))
    catch e
        @error "Failed to read $input" exception=(e, catch_backtrace())
        return (0, 0)
    end
    
    # Check for empty content
    if isempty(nwk)
        @warn "Input file $input is empty, skipping"
        return (0, 0)
    end
    
    tree_repeated_taxa = 0 
    tree_insufficient_taxa = 0 
    
    # Process the single tree
    temp = try
        count_missing(nwk)
    catch e
        @error "Failed to process tree in $input" exception=(e, catch_backtrace())
        return (0, 0)
    end
    
    if isnothing(temp)
        # Tree has repeated taxa - invalid
        println("The input $input has repeated taxa")

        # Remove intermediate mapping files if the tree is invalid 
        # Note: These files may not exist depending on dup_rate/loss_rate settings
        if !debug_mode
            try
                if !isempty(intermediate_file_mapsl) && isfile(intermediate_file_mapsl)
                    rm(intermediate_file_mapsl; force=true)
                end
            catch e
                @warn "Failed to remove $intermediate_file_mapsl" exception=(e, catch_backtrace())
            end

            try
                if !isempty(intermediate_file_maplg) && isfile(intermediate_file_maplg)
                    rm(intermediate_file_maplg; force=true)
                end
            catch e
                @warn "Failed to remove $intermediate_file_maplg" exception=(e, catch_backtrace())
            end

            # Remove the tree file if it is invalid
            try
                if !isempty(input) && isfile(input)
                    rm(input; force=true)
                end
            catch e
                @warn "Failed to remove $input" exception=(e, catch_backtrace())
            end 
        else
            println("  DEBUG: Keeping intermediate files and tree for analysis")
        end

        tree_repeated_taxa = 1
        
    elseif temp[1] > max_taxa_missing 
        # Tree has too many missing taxa - invalid
        println("""The input $input has too many taxa missing (temp[1]=$(temp[1]), max_taxa_missing=$max_taxa_missing)""")

        # Remove intermediate mapping files if the tree is invalid
        if !debug_mode
            try
                if !isempty(intermediate_file_mapsl) && isfile(intermediate_file_mapsl)
                    rm(intermediate_file_mapsl; force=true)
                end
            catch e
                @warn "Failed to remove $intermediate_file_mapsl" exception=(e, catch_backtrace())
            end
            
            try
                if !isempty(intermediate_file_maplg) && isfile(intermediate_file_maplg)
                    rm(intermediate_file_maplg; force=true)
                end
            catch e
                @warn "Failed to remove $intermediate_file_maplg" exception=(e, catch_backtrace())
            end
            
            # Remove the tree file if it is invalid
            try
                if !isempty(input) && isfile(input)
                    rm(input; force=true)
                end
            catch e
                @warn "Failed to remove $input" exception=(e, catch_backtrace())
            end
        else
            println("  DEBUG: Keeping intermediate files and tree for analysis")
            println("  DEBUG: Tree content: $nwk")
            println("  DEBUG: Modified would be: $(temp[2])")
        end 

        tree_insufficient_taxa = 1
        
    else
        # Tree is valid - write the modified tree to output
        modified_tree = temp[2]
        try
            open(output, "w") do io
                write(io, modified_tree * "\n")
            end
        catch e
            @warn "Failed to write output $output" exception=(e, catch_backtrace())
        end
    end
    
    # Return counts of invalid trees
    return tree_repeated_taxa, tree_insufficient_taxa  
end


"""
    modify_newick_for_n_genes(input::String, output::String, n::Int, iteration::Int) 
        -> Tuple{Int, Int}

Run `modify_newicks` across multiple gene tree files for a given replicate and iteration.

This function processes `n` gene tree files, renaming valid trees and writing them 
to the specified output directory. It filters out trees with:
- repeated taxa (non-hidden paralogy),
- fewer than 4 remaining taxa (i.e., more than 4 missing from the set A–H).

# Arguments
- `input::String`: Input directory containing original Newick gene tree files.
- `output::String`: Output directory where modified trees will be written.
- `n::Int`: Number of gene tree files to process.
- `iteration::Int`: Replicate or iteration index (used in output file naming).

# Returns
- `Tuple{Int, Int}`: Total number of excluded trees due to:
    1. Repeated taxa (non-hidden paralogy),
    2. Too many missing taxa.

# Notes
Each gene tree file is expected to be named `g_trees{GENEID}.trees`, where 
`GENEID` is zero-padded to match the total number `n`. Output files are named 
using the pattern:
    `g_trees_noLocusID_Gene{GENEID}_Int{iteration}.trees`.

# Examples
```julia-repl
julia> modify_newick_for_n_genes("input_dir", "output_dir", 100, 1)
(8, 21)  # 8 trees removed for repeated taxa, 21 for missing taxa
```
"""
function modify_newick_for_n_genes(input::String, 
                                output::String, 
                                n::Int, 
                                iteration::Int, 
                                dup_rate::Float64,
                                n_inds::Int, 
                                max_taxa_missing::Int,
                                debug_mode::Bool=false
                                )
    
    # report the number of trees with repeated taxa and insufficient taxa 
    # Those trees are removed from the dataset 
    tree_repeated_taxa_tot = 0 
    tree_insufficient_taxa_tot = 0 

    # Below are statistics that are trees kept for downstream analysis 
    # Report the number of trees experiencing gene loss and gene duplication
    num_trees_experiencing_gene_duplication_and_loss = 0
    # num of trees only experiencing gene loss but no gene duplication
    num_trees_experiencing_gene_loss_only = 0
    # num of trees experiencing nothing (i.e. no gene loss and no gene duplication)
    num_trees_experiencing_nothing = 0

    # count how many leaf left in trees 
    # Later, will be used to calculate the average number of taxa in the final dataset 
    num_leaf_left_in_trees = Int[] 
    
    # Save a list of gene trees name for each category 
    gene_trees_gene_duplication_and_loss = String[]
    gene_trees_gene_loss_only = String[]
    # No need to save for gene tree experiencing nothing
    # because the rest are gene trees experiencing nothing

    for gene_tree in 1:n
        genenum_string = pad_number(gene_tree, n) # see Part 3 "pad_number" below
        output_tree_name = "g_trees_noLocusID_Gene$(genenum_string)_Int$(iteration).trees"
        input_dir = joinpath(input, "g_trees$genenum_string.trees")
        output_dir = joinpath(output, output_tree_name)

        # While simulating, deleting temp files from trees with repeated taxa and insufficient taxa
        # Two files need to be deleted (only if the tree is invalid):
        # We want to only keep the intermediate files for the valid trees 
        # Then we can analyze it to see how many trees experiencing gene loss only 
        # and how many trees experiencing (gene loss + gene duplication) or experiencing nothing  
        intermediate_file_mapsl = joinpath(input, "$genenum_string.mapsl")
        intermediate_file_maplg = joinpath(input, "$(genenum_string)l1g.maplg")

        tree_repeated_taxa, tree_insufficient_taxa = modify_newicks(input_dir, 
                                                                    output_dir, 
                                                                    intermediate_file_mapsl, 
                                                                    intermediate_file_maplg,
                                                                    max_taxa_missing,
                                                                    debug_mode
                                                                    )

        tree_repeated_taxa_tot += tree_repeated_taxa
        tree_insufficient_taxa_tot += tree_insufficient_taxa

        #= Analyze the files after each gene tree modification
        # For the .mapsl file left after filtering: 
        # It could be three situations: 
        # 1) Only gene loss happenes and no gene duplication 
            -> gene loss only but no gene duplication 
            -> No "dup" in .mapsl file and number of "leaf" word < eight 
            -> also count how many taxa are missing
            -> Count how many times loss happens 
        # 2) Gene duplication happens followed by gene loss 
            -> This could be three senarioes (see simulation_postprocess.jl and below) 
            -> The word "dup" in .mapsl file 
            -> also count if and if how many taxa are missing
            -> Count how many times loss and dup happens 
                A few senarios could happen here: 
                    The below cases will be handled during post-processing steps:  
                    a. Gene loss happened after duplication and the gene copy 
                    got duplicated and then lost. In this case, no true paralogy 
                    is present, but the gene copy is lost 
                        --> This is still not a hidden paralogy event
                    In this case, I examine if all remaining gene copies
                    are from the same accession in the locus tree. If so, this is still not a hidden paralogy event.
                    b. Weak hidden paralogy: Gene duplication happened and then got lost. However, 
                    locus tree has the same topology as the species tree, 
                    the hidden paralogy changes no topology but it changes the 
                    branch length 
                    c. Strong hidden paralogy: the remainning gene copies
                    are generated from duplications and the topology is changed. 
        # 3) No gene loss and no gene duplication 
            -> no gene duplication and no gene loss 
            -> Should have eight "leaf"
        # Below code will check the .mapsl file and 
        determine if the above three situations happened =# 
        
        # If dup_rate == 0, no mapsl file so skip the below analysis 
        if dup_rate > 0 
            # When duplicate >= 1 and/or loss >= 1: 
            # Both mapsl and maplg files should exist 
            # When duplicate == 0 and loss >= 0: 
            # only maplg file should exist 
            if !isempty(intermediate_file_mapsl) && isfile(intermediate_file_mapsl) 
                
                # Read lines and change everything into lower case 
                lines = lowercase.(readlines(intermediate_file_mapsl)) 

                # Calculate the number of "leaf", "loss" and "dup" in the .mapsl file 
                num_leaf = count(line -> occursin(r"(?i)leaf", line), lines)
                num_loss = count(line -> occursin(r"(?i)loss", line), lines)
                num_dup = count(line -> occursin(r"(?i)dup", line), lines)

                # situation 1: only gene loss happens 
                if num_dup == 0 && num_loss >= 1 
                    # update the number of trees experiencing gene loss only 
                    num_trees_experiencing_gene_loss_only += 1 
                    push!(num_leaf_left_in_trees, num_leaf) 
                    push!(gene_trees_gene_loss_only, output_dir)

                # situation 2: gene duplication happens 
                # gene loss and gene duplication happened
                elseif num_dup >= 1 && num_loss >= 1
                    num_trees_experiencing_gene_duplication_and_loss += 1
                    push!(num_leaf_left_in_trees, num_leaf) 
                    push!(gene_trees_gene_duplication_and_loss, output_dir)

                # situation 3: no gene loss and no gene duplication 
                # meaning no gene duplication and no gene loss 
                elseif num_dup == 0 && num_loss == 0
                    num_trees_experiencing_nothing += 1 
                    push!(num_leaf_left_in_trees, num_leaf) 

                else
                    error("Unexpected case in $intermediate_file_mapsl check!")
                end
            end 
        else # if dup_rate == 0 
            # Gene trees will never experience gene loss and gene duplication
            num_trees_experiencing_nothing += 1 
            push!(num_leaf_left_in_trees, 8 * n_inds)
        end 

    end


    # This handles nan as well 
    avg_num_leaf_left_in_trees = isempty(num_leaf_left_in_trees) ? NaN : mean(num_leaf_left_in_trees) 
    println("Average number of leaf left in trees (one rep): $avg_num_leaf_left_in_trees")

    #------ Finally, clean up intermediate files -------# 
    # remove intermediate files to save space under non-debugging mode  
    if debug_mode == false # debugging mode will keep intermediate files 
        
        # check if any files ending with ".mapsl" or ".maplg" exist in the input directory 
        # Check if any .mapsl or .maplg files exist in the input directory.
        # If any are present, keep this folder (no-op); otherwise the existing
        # `else` block will remove the parent iteration folder.
        if !isempty(collect(glob("*.maplg", input))) # ONLY search maplg since dup=0 means only maplg exists 
            # There are intermediate mapping files to keep
            # -> only remove the intermediate files that are not needed 
            # Under debug mode, only save g_trees and l_trees files 
            files_to_remove = [".mapsl", ".maplg", "_tree.trees"] 
            for suffix in files_to_remove
                pattern = "*$suffix" 
                for file in glob(pattern, input)
                    rm(file)
                end
            end

            filepath_to_remove = joinpath(input, "../") 
            file_suffix_to_remove = ["command", "conf", "db", "params"] 
            for suffix in file_suffix_to_remove
                pattern = "*.$suffix"
                for file in glob(pattern, filepath_to_remove)
                    rm(file)
                end
            end

            file_prefix_to_remove = "simphysim"
            pattern = "$file_prefix_to_remove*"
            for file in glob(pattern, filepath_to_remove)
                rm(file)
            end

        else 
            # if there is no mapsl or maplg file, we can remove this interation folder 
            # remove the previous folder layer of input directory 
            parent_dir = dirname(input) # This is the interation folder 
            if isdir(parent_dir) 
                rm(parent_dir; force=true, recursive=true)
            end 
        end 

        # Move output files one level up to remove nested Int/1 folder structure
        # Find files in the output_dir that end with .trees and move them up one level
        if isdir(input)
            # Get the parent directory (iteration folder)
            iteration_dir = dirname(input)
            
            # Find all .trees files in output_dir
            trees_files = filter(f -> endswith(f, ".trees"), readdir(input))
            
            for file in trees_files
                src_path = joinpath(input, file)
                dst_path = joinpath(iteration_dir, file)
                # Move file from nested location to iteration folder
                mv(src_path, dst_path; force=true)
            end
            
            # Remove the now-empty nested folder if it's empty
            if isempty(readdir(input))
                rm(input; recursive=true)
            end
        end
    end 

    #= return: 
    1. number of trees with repeated taxa 
    2. number of trees with insufficient taxa 
    3. number of trees experiencing gene loss only 
    4. number of trees experiencing gene loss and gene duplication
        -> Later, during post-processing steps, will check if the remaining
              gene copies are from the same accession or not
              to determine if it is a true hidden paralogy event 
    5. number of trees experiencing nothing
    6. average number of leaf left in the final tree dataset 
    =# 
    return tree_repeated_taxa_tot, 
        tree_insufficient_taxa_tot, 
        num_trees_experiencing_gene_loss_only, 
        num_trees_experiencing_gene_duplication_and_loss, 
        num_trees_experiencing_nothing, 
        avg_num_leaf_left_in_trees,
        gene_trees_gene_loss_only,
        gene_trees_gene_duplication_and_loss
end 

#-----------------------------------------------#       
#               Part 2 
#    Re-run simphy until:
# 1) n_genes hit the target; or 
# 2) max_iteration reaches   
#-----------------------------------------------#
"""
    seed_generator(master_seed::Int, n::Int, m::Int, 
                    output_dir::String, output_file_name::String) 
        -> Matrix{Int32}

Generate a matrix of integer seeds for SimPhy simulations and write them to a file.

Each seed is drawn from a `StableRNG` initialized with `master_seed`.
A total of `n × m` seeds are generated, where `n` is the number of 
replicates and `m` is the number of iterations per replicate.

# Arguments
- `master_seed::Int`: Master seed to start the random number generator (`Xoshiro`).
- `n::Int`: Number of replicates (`n_reps`).
- `m::Int`: Number of iterations per replicate (`max_iteration`).
- `output_dir::String`: Directory where the output file will be written.
- `output_file_name::String`: File name to write the seeds (e.g. `"random_seeds.txt"`).

# Returns
- `Matrix{Int32}`: A matrix of size `(n, m)` containing the generated seeds.

# Output File Format
- First line: column headers (`col1`, `col2`, ..., `colm`)
- Subsequent lines: one line per rep, with `m` seed values per line (tab-separated)

# Examples
```julia-repl
julia> seeds = seed_generator(42, 3, 4, ".", "random_seeds.txt")
3×4 Matrix{Int32}:
 1313536569  1702490664  1408128490  1063381724
 1013305710  1694460136  1777051044   777323194
  820768661  1426708740  1792281826  1724743702
```
"""
function seed_generator(master_seed::Int, n::Int, m::Int, 
                        output_dir::String, output_file_name::String) 

    rng = StableRNG(master_seed)
    # generates n x m seed array with positive numbers only: 
    seeds = abs.(rand(rng, Int32, n, m)) 

    # write all seeds to a random_seeds.txt file in the output_dir 
    output_file = joinpath(output_dir, output_file_name)
    open(output_file, "w") do io
        for j in 1:m
            print(io, j == m ? "col$j\n" : "col$j\t") 
            # Above: if j is the last then separate by a new line
        end

        for i in 1:n
            for j in 1:m
                print(io, j == m ? "$(seeds[i, j])\n" : "$(seeds[i, j])\t") 
                # Above: if j is the last then separate by a new line
            end
        end
    end
    return seeds
end  


"""
    calculate_batch(num::Int, n_genes::Int) -> Int

Estimate the number of gene trees (batch size) to simulate to reach the target count, 
given the number of existing gene trees.

This function assumes a missing rate based on how many of the `n_genes` are present.
It then estimates how many additional trees need to be simulated to compensate 
for expected missing ones, rounding up to the nearest integer.
However, the estimated batch size will not exceed 3 times of the n_genes
since otherwise it may lead to excessive resource usage.

# Arguments
- `num::Int`: Number of gene trees currently present.
- `n_genes::Int`: Total number of target gene trees.

# Returns
- `Int`: Estimated batch size to simulate next.

# Calculation
Let `missing = n_genes - num` and `success_rate = num / n_genes`.  
Then: `estimated_batch = ceil(missing / success_rate)`.  
If `success_rate == 0`, return `n_genes`.

# Examples
```julia-repl
julia> calculate_batch(24, 30)  # 20% missing rate
8

julia> calculate_batch(0, 30)   # All genes missing
30
```
"""
function calculate_batch(num::Int, n_genes::Int) 
    num_missing = n_genes - num 
    success_rate = num / n_genes 
    if success_rate == 0
        estimated_batch = n_genes 
    elseif success_rate <= 0.2
        println("success rate too small, estimated batch size = 5 * n_genes")
        estimated_batch = n_genes * 6
        # This limits the number of genes got simulated each time 
        # Or if success_rate is too small, the estimated batch size will be too big 
        # This will explode the storage
    else 
        estimated_batch = round(num_missing / success_rate) 
    end 
    # println("For replicate", n_reps, " estimated batch size = ", estimated_batch)
    return Int(estimated_batch)
end 


#-----------------------------------------------#       
#               Part 3  
#       Other utilities functions     
#-----------------------------------------------#
"""
    pad_number(num::Int, range::Int) -> String

Convert `num` into a zero-padded string, using the number of digits required to represent `range`.

This is useful for naming consistency, e.g. when generating file names for SimPhy gene trees
like `"g_trees001.trees"`, `"g_trees002.trees"`, etc.

# Arguments
- `num::Int`: The number to pad.
- `range::Int`: The maximum possible value in the range, which determines the padding width.

# Returns
- `String`: A zero-padded string version of `num`.

# Examples
```julia-repl
julia> pad_number(1, 10)
"01"

julia> pad_number(1, 5)
"1"

julia> pad_number(1, 200)
"001"
```
"""
function pad_number(num::Int, range::Int) 
  digits = ceil(Int, log10(range + 1))
  number_string = lpad(string(num), digits, '0') 
  return number_string
end

"""
    check_existing_dir(dir_input::Vector{String})

Interactively check if one or more directories exist, and optionally remove them.

This function is useful for testing simulation code repeatedly by letting the user
decide whether to remove existing output folders before proceeding.

# Behavior
- If the path(s) exist:
  - For multiple paths: prompts the user to remove all (`ALL`) or remove individually (`N`).
  - For individual paths: prompts the user to remove (`Y`) or not (`N`).
- If the user chooses not to remove a path, an error is thrown to halt execution.
- If a path does not exist, it is ignored.
- If the input vector is empty, an error is raised.

# Arguments
- `dir_input::Vector{String}`: A list of one or more directory paths to check.

# User Inputs
- `ALL` / `all`: Remove all directories (if multiple are provided).
- `Y` / `y`: Remove an individual directory.
- `N` / `n`: Skip removal of a directory (raises an error).
- Any other input: Will be reprompted.

# Examples
```julia-repl
julia> check_existing_dir(["sim_output1", "sim_output2"])
Multiple directories were provided:
Do you want to remove all directories at once? ALL/N
```
"""
function check_existing_dir(dir_input::Vector)

    len = length(dir_input)
    if len == 0 # If input an empty string, then output an error 
        error("Oh no! It's empty string/vector. Check out the code.")
        return 
    end

    if len > 1 && ispath(dir_input[1]) 
        # if input is a list of strings and the first path of the list exists 
        println("Multiple directories were provided: ")
        while true
            println("Do you want to remove all directories at once? ALL/N")
            user_input = lowercase(readline())
            if user_input == "all"
                for dir in dir_input
                    if ispath(dir)
                        rm(dir; recursive=true) # Remove all paths at once 
                    else
                        println("Directory $dir does not exist.")
                        continue
                    end 
                end
                println("All directories have been removed.") 
                return    
            elseif user_input == "n"
                println("Don't want to remove all dirs, let's remove each individual one")
                break 
            else
                println("Invalid input. Please type ALL or N.")
            end
        end
    end

    for dir in dir_input # remove individual strings in the input
        # If only one path in dir_input -> 
        if ispath(dir)
            while true
                println("The directory $dir already exists. Do you want to remove it? Y/N")
                user_input = lowercase(readline()) 
                if user_input == "y"
                    rm(dir; recursive=true)
                    break
                elseif user_input == "n"
                    error("Hey! $dir already exists, and we don't want to remove it")
                    return
                else
                    println("Invalid input. Please type Y or N.")
                end
            end
        end
    end
end

"""
    replace_tips_with_letters(tree::Union{HybridNetwork, String}) 
        -> Union{HybridNetwork, String}

Replace tip labels in a tree with capital letters A–Z.

This function simplifies tip names in either a Newick string or a 
`HybridNetwork` tree object. It supports up to 26 taxa, assigning letters 
`A`, `B`, `C`, ..., based on their order of appearance. If more than 26 
tips are present, the original input is returned with a warning.

# Arguments
- `tree::Union{HybridNetwork, String}`: Either a Newick-formatted string or 
  a `HybridNetwork` object.

# Returns
- `HybridNetwork`: Modified tree with simplified tip labels if input was a 
  tree.
- `String`: Modified Newick string with simplified tip labels if input was 
  a string.

# Notes
- Tip names are replaced based on order of appearance, not taxonomic 
  identity.
- The function handles special characters like `*` and branch lengths 
  correctly.

# Examples
```julia-repl
julia> replace_tips_with_letters(
           "(Homo:0.01,(((Cat:0.004,Dog:0.003):0.007,(Dinosaur:0.02,"
           * "Birds:0.02):0.01)));")
"(A:0.01,(((B:0.004,C:0.003):0.007,(D:0.02,E:0.02):0.01)));"

julia> tree = readnewick(
           "(Homo:0.01,(((Cat:0.004,Dog:0.003):0.007,(Dinosaur:0.02,"
           * "Birds:0.02):0.01)));")
julia> replace_tips_with_letters(tree)
HybridNetwork, Rooted Network
9 edges, 10 nodes (5 tips, 5 internal).
Tip labels: A, B, C, D, ...

julia> replace_tips_with_letters(
           "(Homo:0.01*4,(((Cat:0.004*5,Dog:0.003*10):0.007,"
           * "(Dinosaur:0.02*20,Birds:0.02*3):0.01*9)));")
"(A:0.01*4,(((B:0.004*5,C:0.003*10):0.007,(D:0.02*20,E:0.02*3):0.01*9)));"
```
"""
function replace_tips_with_letters(tree::Union{HybridNetwork, String}) 
	#= Goal: replace species tips with A, B, C, etc, 
        (but tips should be less than or equal to 26). 
    The input (tree) can both be a string or a tree
    If input is tree, then return a tree object
    If input is a string, then return a string: 
        which could be used to modified the input to SimPhy=#

    new_labels = ['A':'Z';]
    if isa(tree, HybridNetwork) # If input is a tree, output a tree 
        tips = tiplabels(tree)

        if length(tips) > 26 
            #Give an warning and return the original tree if there are more than 26 tips 
            println("Warning: The number of tips exceeds 26. That's too much. Exiting.")
            return tree
        end 
    
        for (i, tip) in enumerate(tree.leaf)
            tip.name = string(new_labels[i])
        end
        return tree
    
    elseif isa(tree, String) # if input is a newick string, output a string
        tips = split(tree, [',', '(', ')', ':', '*', ';'])  # Extract tips
        tips = filter(x -> !(x in ["", ":", ";"]) && !occursin(r"\d", x), tips) 
        # Above -> Remove empty, :, and any string containing numbers

        if length(tips) > 26
            println("Warning: The number of tips exceeds 26. That's too much. Exiting.")
            return tree
        end
        for (i, tip) in enumerate(tips)
            tree = replace(tree, tip => string(new_labels[i]))
        end
        return tree
    else # Give a warning if the input is neither network nor string 
        println("Warning: The input needs to be a HybridNetwork or String.")
        return tree
    end
end


"""
    setup_rep_output_folders(
        folder_path_list::Vector{String}, 
        simulation_rep::Int, 
        path_string::String
    ) -> String

Construct an output path for a simulation replicate by appending a string to 
the corresponding folder path.

# Arguments
- `folder_path_list::Vector{String}`: A list of folder paths, one per replicate.
- `simulation_rep::Int`: The current replicate ID (1-based index).
- `path_string::String`: The subdirectory or filename to append.

# Returns
- `String`: The constructed output path using `joinpath`.

# Examples
```julia-repl
julia> folder_list = ["rep1_dir", "rep2_dir", "rep3_dir"];

julia> setup_rep_output_folders(folder_list, 2, "results.txt")
"rep2_dir/results.txt"
```
"""
function setup_rep_output_folders(
    folder_path_list::Vector, 
    simulation_rep::Int, 
    path_string::String
    )
  rep_folder_path = folder_path_list[simulation_rep]  
  return (joinpath(rep_folder_path, path_string))
end 

"""
    write_boostrapTrees_2_Ind_filepath(
        file_path::String, 
        output_dir::String, 
        output_filename::String, 
        relative_path::Bool
    )

Split a tree file with multiple bootstrap trees into individual files and write 
their paths to a text file.

This function is used for testing. It is not part of the main pipeline but can 
help prepare bootstrap trees in a format compatible with bootsnaq, where input 
should be formatted as a vector of individual tree file paths.

# Arguments
- `file_path::String`: Path to a file containing multiple Newick bootstrap trees 
  (one per line).
- `output_dir::String`: Directory where individual tree files will be written.
- `output_filename::String`: Name of the output text file that will store the 
  list of individual tree file paths.
- `relative_path::Bool`: Whether to write tree paths as relative (`true`) or 
  absolute (`false`).

# Behavior
- Each line in `file_path` is treated as a single tree and written to a new file.
- File paths to all individual tree files are written to `output_filename`.

# Returns
- Nothing is returned. Files are written to disk.

# Examples
```julia-repl
julia> write_boostrapTrees_2_Ind_filepath(
           "bootstrap_trees.txt", 
           "trees_split/", 
           "tree_list.txt", 
           true
       )
```
"""
function write_boostrapTrees_2_Ind_filepath(
        file_path::String, 
        output_dir::String, 
        output_filename::String, 
        relative_path::Bool
        )
    trees = readlines(file_path)
    output_txt = joinpath(output_dir, output_filename)

    open(output_txt, "w") do io 
        for (i, tree) in enumerate(trees) 
            if relative_path
                tree_path = joinpath(output_dir, "tree_$i.trees")
            else 
                tree_path = joinpath("./", "tree_$i.trees") 
            end 
            open(tree_path, "w") do f
                write(f, tree)
            end 
            println(io, tree_path) 
        end 
    end
end

"""
    replace_pop_inInd_file(input_file::String, output_file::String)

Replace `"POP"` in the third column of a `.ind` file with a group label extracted 
from the first column.

This function is used to process `.ind` files where "POP" is a placeholder for 
population labels. It assigns a population letter based on the first part of the 
ID in column 1 (e.g., `"A_0"` or `"A_1"` → `"A"`).

# Arguments
- `input_file::String`: Path to the input `.ind` file. Each line must contain 
  three tab-separated fields.
- `output_file::String`: Path to the output file where the updated lines will be 
  written.

# Returns
- Nothing. A new `.ind` file is written to `output_file` with updated population 
  labels.

# Input Format
Each line must have exactly three tab-separated columns:
1. Accession ID (e.g., `A_0`, `B_1`)
2. Numeric value
3. "POP" (which will be replaced)

# Example
Input line:
"""
function replace_pop_inInd_file(input_file::String, output_file::String)
    open(input_file, "r") do infile
        open(output_file, "w") do outfile
            for line in eachline(infile)
                fields = split(line, '\t')
                if length(fields) != 3
                    error("Check $(input_file) -> .ind format should have three columns")
                end 
                if fields[3] != "POP"
                    error("Check $(input_file) -> .ind format has last column in \'POP\'")
                end 
                fields[3] = split(fields[1], '_')[1]  # Extract the letter before "_"
                println(outfile, join(fields, '\t'))
            end
        end
    end
end

"""
    generate_software_seeds(
        master_seed::Int, 
        software_names::Vector{String}
    ) -> Dict{String, Int}

Generate a dictionary of stable random seeds for each software name using a 
master seed.

This function produces reproducible seeds for software components in a pipeline, 
based on a given `master_seed`. The same input will always result in the same 
set of seeds due to the use of `StableRNG`.

# Arguments
- `master_seed::Int`: The initial seed used to initialize the random number 
  generator.
- `software_names::Vector{String}`: A list of software names that each require 
  a unique seed.

# Returns
- `Dict{String, Int}`: A dictionary mapping each software name to a unique 
  stable integer seed.

# Examples
```julia-repl
julia> master_seed = 123456;

julia> software_names = ["SimPhy", "RAxML", "Astral", "SeqGen", "snaq"];

julia> generate_software_seeds(master_seed, software_names)
Dict("SimPhy" => 1806341205, 
    "RAxML" => 1735483215, 
    "Astral" => 1170911113, 
    "SeqGen" => 1313482870, 
    "snaq" => 1234567890)
"""
function generate_software_seeds(master_seed::Int, software_names::Vector{String})
    rng = StableRNG(master_seed) 
    seeds = abs.(rand(rng, Int, length(software_names)))  
    return Dict(software_names[i] => seeds[i] for i in eachindex(software_names))
end

"""
    generate_master_seed(params::Dict{String, Any}) -> Int
Generate a deterministic master seed based on parameter values by encoding 
them into products of prime numbers.

# Description
- **Continuous values** (e.g., `0.01`, `0.0001`) are encoded by assigning digits 
  after the decimal point to specific prime numbers.
- **Categorical values** (e.g., substitution models like `"G"` or `"LG"`) are 
  mapped to fixed primes using a hardcoded lookup (`category_map`).
- **Integer values** are mapped to the next largest prime using `nextprime()`.

The final seed is the product of all prime values assigned to each parameter.

# Arguments
- `params::Dict{String, Any}`: A dictionary where keys are parameter names and 
  values are numeric or categorical values.

# Returns
- `Int`: A stable integer seed representing the combined encoded values of the 
  parameters.

# Notes
- The `digit_to_prime` map converts each digit character to a unique prime.
- Category mappings (e.g., `"ratevar"` → `"G"` or `"LG"`) are hardcoded inside 
  the function.
- Internally uses `big()` to avoid overflow, then converts to `Int`.

# Examples
```julia-repl
julia> params = Dict("dup" => 0.01, "loss" => 0.01, "ratevar" => "G");

julia> generate_master_seed(params)
14613783
```
"""
function generate_master_seed(params::Dict{String, Any})
    # println(params)
    digit_to_prime = Dict(
        '0'=>2, '1'=>3, '2'=>5, '3'=>7, 
        '4'=> 11, '5'=>13, '6'=>17, '7'=>19, 
        '8'=>23, '9'=>29, 
        '.' => 31, 'e' => 37, '-' => 41)
    category_map = Dict("ratevar" => Dict("N"=>43, "G"=>47, "L" => 59, 
                                        "LG"=>53, "GL" => 53, "G*L" => 53)) 
    master_seed = big(1)  # big to prevent large integer overflow 
    prime = 1 # prime starts from 1 

    for (param, value) in params
        
        # Skip parameters with nothing value (e.g., genelen when not provided)
        if value === nothing
            continue
        end
        
        if haskey(category_map, param) # If parameter in category_map 
            prime = category_map[param][value]

        elseif isa(value, Float64) || isa(value, Float32) # parameter = number  
            splitted_number = split(string(value), '.')

            if length(splitted_number) > 1 
                digits_after_decimal_place = splitted_number[2]  
                for ch in digits_after_decimal_place
                    prime *= digit_to_prime[ch]
                end
            else 
                digits_before_decimal_place = splitted_number[1]
                for ch in digits_before_decimal_place
                    prime *= digit_to_prime[ch] 
                end 
            end

        elseif isa(value, Int)
            prime = nextprime(value) 
            # Above: nextprime assigns a prime number which is bigger than the value
        else
            error("""Oh no something is wrong with the parameters setting. 
                    Unable to generate master seed""")
        end   
        master_seed *= prime
    end
    return Int(master_seed) 
    # Above: change to integer because generate_software_seeds only takes integers 
end

"""
    sample_substitution_params(seed)
Return a stable set of parameters for a given seed.
Input: 
- `seed`: Seed for the random number generator 
Output: 
- `kappa::Float64`: Transition/transversion ratio sampled from LogNormal(1.4215, 0.2798).
- `basefreqs::Vector{Float64}`: Base frequencies sampled from Dirichlet(66.59, 38.41, 38.61, 67.12).
- `alpha::Float64`: Gamma shape parameter sampled from Gamma(3.267, 0.109).
# Examples
```julia-repl
julia> sample_substitution_params(42)
(kappa = 3.4348167506335687, basefreqs = [0.31818654496208937, 0.21258062786139054, 0.19404110221414664, 0.27519172496237343], alpha = 0.628276797912641)
``` 
This information is based on the readme file: 
- to simulate each gene with its own substitution model, use HKY with:
  * kappa from LogNormal(μ=1.4215, σ=0.2798)
  * frequencies from Dirichlet(66.59, 38.41, 38.61, 67.12)
  * alpha from Gamma(α=3.267, θ=0.109). 
"""
function sample_substitution_params(seed) 
    rng = StableRNG(seed)   # stable RNG object
    # κ ~ LogNormal
    kappa = rand(rng, LogNormal(1.4215, 0.2798))
    # base frequencies ~ Dirichlet
    freq_dist = Dirichlet([66.59, 38.41, 38.61, 67.12])
    basefreqs = rand(rng, freq_dist)
    # γ-shape α ~ Gamma
    alpha = rand(rng, Gamma(3.267, 0.109))
    return (kappa=kappa, basefreqs=basefreqs, alpha=alpha)
end

"""
    get_dict_for_seed_setting(params_string::String) -> Dict{String, Any}

Parse a simulation parameter string into a dictionary used for master seed 
generation.

This function extracts values for duplication rate, loss rate, rate variation 
model, and number of individuals from a specially formatted string. The 
dictionary it returns can be passed directly to `generate_master_seed`.

# Expected Format
The input string should be in the format: 
    DUP<dup_rate>-LOS<loss_rate>-RV<ratevar>-N_ind<n_inds>  

Examples:
- `"DUP0.01-LOS0.02-RVG-N_ind10-genelen1000"`
- `"DUP1.0e-5-LOS0.0-RVLG-N_ind50-genelen500"`

# Arguments
- `params_string::String`: A parameter string with values for `dup_rate`, 
  `loss_rate`, `ratevar`, and `n_inds`, joined by labeled segments.

# Returns
- `Dict{String, Any}`: A dictionary with the following keys:
    - `"dup_rate"`: `Float32` — duplication rate
    - `"LOSS_RATE"`: `Float32` — loss rate
    - `"ratevar"`: `String` — rate variation model (e.g. `"G"`, `"LG"`)
    - `"n_inds"`: `Int` — number of individuals
    - `"scaling_factor"`: `Float32` — scaling factor for seed generation
    - `"genelen"`: `Float32` or `nothing` — gene length if provided, else `nothing` 

# Notes
- Parsing is done using string splitting at fixed substrings: `-LOS`, `-RV`, 
  and `-N_ind`.
- This approach avoids regex due to issues with scientific notation (e.g., 
  `1.0e-5`).

# Examples
```julia-repl
julia> get_dict_for_seed_setting("DUP0.01-LOS0.02-RVG-N_ind10-genelen1000")
Dict("dup_rate" => 0.01f0,
     "LOSS_RATE" => 0.02f0,
     "ratevar" => "G",
     "n_inds" => 10,
     "scaling_factor" => 1000,)
```
""" 
function get_dict_for_seed_setting(params_string::String)

    # Split based on known breakpoints
    part1 = split(params_string, "-LOS")[1] # contains DUP rate
    rest1 = split(params_string, "-LOS")[2]

    part2 = split(rest1, "-RV")[1] # contains LOSS rate
    rest2 = split(rest1, "-RV")[2]

    part3 = split(rest2, "-N_ind")[1] # contains RV
    rest4 = split(rest2, "-N_ind")[2] # contains n_inds + SF

    part4 = split(rest4, "-SF")[1] # contains n_inds 
    part5 = split(rest4, "-SF")[2] # contains SF + optional genelen
    
    # Some legacy trails might not contain genelen info 
    has_genelen = occursin("-genelen", part5)
    if has_genelen
        scaling_factor_str = split(part5, "-genelen")[1] # extract SF before genelen
        part6 = split(part5, "-genelen")[2] # contains genelen 
    else
        scaling_factor_str = part5
    end

    dup_rate_str = replace(part1, "DUP" => "")
    loss_rate_str = part2
    ratevar = part3
    n_inds_str = part4

    # Convert to float and int 
    dup_rate = parse(Float32, dup_rate_str)
    loss_rate = parse(Float32, loss_rate_str)
    n_inds = parse(Int, n_inds_str)
    scaling_factor = parse(Float32, scaling_factor_str) 

    if has_genelen
        genelen = parse(Float32, part6)
    else
        genelen = nothing
    end

    return Dict(
        "dup_rate" => dup_rate,
        "LOSS_RATE" => loss_rate,
        "ratevar" => ratevar,
        "n_inds" => n_inds, 
        "scaling_factor" => scaling_factor,
        "genelen" => genelen 
    )
end


"""
    median_branch_length(tree::HybridNetwork) -> Float64

Compute the median branch length of a phylogenetic tree or network.

# Arguments
- `tree::HybridNetwork`: A phylogenetic tree or network object from the 
  PhyloNetworks package.

# Returns
- `Float64`: The median of all non-missing branch lengths in the tree.

# Examples
```julia-repl
julia> using PhyloNetworks
julia> tree = readnewick("(A:0.1,(B:0.2,C:0.3):0.4);");
julia> median_branch_length(tree)
0.25
```
"""
function median_branch_length(tree::HybridNetwork)
    branch_lengths = [e.length for e in tree.edge if !ismissing(e.length)]
    return median(branch_lengths)
end

"""
    transform_newick_no_lineage_tree(
        tree::String, 
        Ne::Int, 
        rate_str::String
    ) -> String

Transform a Newick tree string by scaling branch lengths based on a constant 
effective population size (Ne) and appending a mutation rate string.

This function is used in testing and is not part of the main simulation pipeline.

Each branch length `t` in the Newick string is replaced with `round(2Ne * t)`, 
followed by the given `rate_str` (e.g. `*0.01`). The pattern `:([0-9.]+)` is 
used to identify all branch lengths in the input string.

# Arguments
- `tree::String`: A Newick-formatted tree string.
- `Ne::Int`: Effective population size, used to scale branch lengths.
- `rate_str::String`: String to append after each scaled branch length, 
  typically representing a per-generation mutation rate (e.g. `"0.01"`).

# Returns
- `String`: A Newick string with all branch lengths scaled and annotated.

# Examples
```julia-repl
julia> tree = "(A:0.01,B:0.02);"

julia> transform_newick_no_lineage_tree(tree, 1000, "0.01")
"(A:20*0.01,B:40*0.01);"
```
"""
function transform_newick_no_lineage_tree(tree::String, Ne::Int, rate_str::String)

    pattern = r":([0-9.]+)" # Match branch lengths in the Newick string
    result = ""
    last_end = 0
    factor = 2 * Ne

    for m in eachmatch(pattern, tree)
        start_idx = m.offset
        end_idx = m.offset + length(m.match) - 1
        cu = parse(Float64, m.captures[1])
        scaled_int = Int(round(cu * factor))
        transformed = ":" * string(scaled_int) * "*" * rate_str
        result *= tree[last_end+1:start_idx-1] * transformed
        last_end = end_idx
    end

    result *= tree[last_end+1:end]
    return result
end

"""
    multiple_tree_length_with_2Ne(tree::String, TwoNe::Int) -> String

Multiply all branch lengths in a Newick tree string by a constant factor (2Ne).

This function scales every branch length in the input Newick string by the 
provided `TwoNe` value (typically `2 * Ne` in coalescent models). The function 
matches branch lengths using the pattern `:[0-9.]+` and replaces each with the 
scaled value.

# Arguments
- `tree::String`: A Newick-formatted tree string with numeric branch lengths.
- `TwoNe::Int`: The scaling factor (usually `2 * Ne`) to multiply all branch 
  lengths by.

# Returns
- `String`: A Newick string with all branch lengths scaled by `TwoNe`.

# Examples
```julia-repl
julia> tree = "(A:0.01,B:0.02);";

julia> multiple_tree_length_with_2Ne(tree, 1000)
"(A:10.0,B:20.0);"
```
"""
function multiple_tree_length_with_2Ne(tree::String, TwoNe::Int64)
    pattern = r":[0-9.]+"
    return replace(tree, pattern => s -> begin
        value = parse(Float64, s[2:end])  # Skip the ":"
        ":" * string(value * TwoNe)
    end)
end

"""
    transform_newick_with_lineage_tree(
        tree_CU_str::String, 
        tree_sub_str::String, 
        Ne::Int
    ) -> String

Transform a Newick tree string by scaling branch lengths using data from a 
coalescent-unit (CU) tree and a substitution-unit tree.

Each branch is rewritten as `:generation*rate`, where:
- `generation = CU_length * 2 * Ne`
- `rate = substitution_length / generation`

This format is useful for coalescent simulations with branch lengths in 
generations and mutation rates per generation.

# Arguments
- `tree_CU_str::String`: A Newick-formatted string with branch lengths in 
  coalescent units.
- `tree_sub_str::String`: A Newick-formatted string with branch lengths in 
  substitution units.
- `Ne::Int`: The effective population size (used in `2 * Ne` scaling).

# Returns
- `String`: A Newick string with each branch length transformed to the format 
  `:generation*rate`.

# Notes
- Both input trees must have the same number of branch lengths.
- Uses a regular expression to extract branch lengths, including scientific 
  notation (e.g., `1.0e-5`).

# Examples
```julia-repl
julia> tree_CU_str = "(A:0.01,B:0.02);";

julia> tree_sub_str = "(A:0.005,B:0.01);";

julia> Ne = 1000;

julia> transform_newick_with_lineage_tree(tree_CU_str, tree_sub_str, Ne)
"(A:20000*0.00025000,B:40000*0.00025000);"
```
"""
function transform_newick_with_lineage_tree(
            tree_CU_str::String, 
            tree_sub_str::String, 
            Ne::Int
            )
    pattern = r":([0-9\.eE+-]+)"

    cu_matches = collect(eachmatch(pattern, tree_CU_str))
    sub_matches = collect(eachmatch(pattern, tree_sub_str))

    if length(cu_matches) != length(sub_matches)
        error("Mismatch in number of branch lengths between CU and substitution trees.")
    end

    result = ""
    last_end = 0

    for i in 1:length(cu_matches)
        m_cu = cu_matches[i]
        m_sub = sub_matches[i]

        start_idx = m_cu.offset 
        # Above: offset finds the start index of the match in the string 
        end_idx = m_cu.offset + length(m_cu.match) - 1
        # println("Match $i from $start_idx to $end_idx")

        cu = parse(Float64, m_cu.captures[1])
        sub = parse(Float64, m_sub.captures[1])

        generation = cu * 2 * Ne
        rate = sub / generation

        transformed = @sprintf(":%d*%.8f", generation, rate)

        result *= tree_CU_str[last_end+1:start_idx-1] * transformed
        last_end = end_idx
    end

    result *= tree_CU_str[last_end+1:end]
    return result
end

"""
    add_outgroup_to_newick(
        original_newick::String, 
        out_branch::Float64
    ) -> String

Add an outgroup branch to a Newick tree string.

This function is primarily used for testing, and ensures that SimPhy can 
process a species tree where an early-diverging outgroup is added. The 
outgroup's branch length is computed as the sum of `out_branch` and the 
first internal branch length in the original tree.

# Arguments
- `original_newick::String`: A Newick-formatted string representing the tree.
- `out_branch::Float64`: The branch length to add between the original tree 
  and the new outgroup root.

# Returns
- `String`: A modified Newick string with a new outgroup `"O"` added as the 
  sister to the original tree.

# Notes
- The new outgroup `"O"` receives a branch length equal to the sum of the 
  first original branch length and `out_branch`.
- The original tree is nested and attached with branch length `out_branch`.

# Examples
```julia-repl
julia> tree = "(A:0.01,B:0.02);";

julia> add_outgroup_to_newick(tree, 0.1)
"(O:0.11,(A:0.01,B:0.02):0.1);"
```
"""
function add_outgroup_to_newick(original_newick::String, out_branch::Float64)
    
    branch_lengths = collect(eachmatch(r"\d+\.\d+", original_newick))

    first_bl = parse(Float64, branch_lengths[1].match)

    outgroup_bl = first_bl + out_branch # This is the outgroup branch length 
    trimmed_newick = endswith(original_newick, ";") ? 
        original_newick[1:end-1] : original_newick
    new_tree = "(" *
        "O:$(round(outgroup_bl, digits=6)),$trimmed_newick:" *
        "$(round(out_branch, digits=6))" * ");"

    return new_tree
end

"""
    prune_tip_from_newick(newick_str::String, prefix::String) -> String

Remove all tips from a Newick tree string whose names start with a given prefix.

This function parses the Newick string, finds all tips whose names begin with 
the provided prefix, prunes them from the tree, and writes the result back to 
Newick format.

# Arguments
- `newick_str::String`: The input Newick-formatted tree string.
- `prefix::String`: The prefix used to identify tip names to be removed.

# Returns
- `String`: A Newick-formatted string with the specified tips removed.

# Notes
- The tree is re-rooted at the first matching tip before pruning.
- Uses `readnewick` and `writenewick` from PhyloNetworks.

# Examples
```julia-repl
julia> tree = "(A1:0.1,B1:0.2,C2:0.3);";

julia> prune_tip_from_newick(tree, "A")
"(B1:0.2,C2:0.3);"
```
"""
function prune_tip_from_newick(newick_str::String, prefix::String)
    tree = readnewick(newick_str)
    tips_to_prune = [leaf.name for leaf in tree.leaf if startswith(leaf.name, prefix)]
    rootatnode!(tree, tips_to_prune[1]) # Ensure the tree is rooted before pruning 
    for tip in tips_to_prune
        deleteleaf!(tree, tip)
    end
    return writenewick(tree) # writenewick replaces the decrepated writeTopology 
end


"""
    map_accessions_to_species_dict(
        gene_trees::Vector{HybridNetwork}, 
        output_dir::String, 
        mode::String = "snaq"
    ) -> Tuple{Dict{String,String}, String}

Generate a mapping from individual accession names to species names, based on 
tip labels in gene trees, and write the mapping file in either SNaQ or ASTRAL 
format.

# Arguments
- `gene_trees::Vector{HybridNetwork}`: A list of gene trees where each tip is 
  labeled as `"Species_Accession"` (e.g., `"A_1"`, `"B_2"`).
- `output_dir::String`: Directory where the output mapping file will be saved.
- `mode::String = "snaq"`: Output file format. Options are:
  - `"snaq"`: CSV file with columns `species` and `individual` 
    (`id_to_species.csv`), for use with `snaq!`.
  - `"astral"`: Space-delimited text file with lines like 
    `"individual species"` (`astral_mapping.txt`), for use with ASTRAL’s `-a` flag.

# Returns
- `Tuple{Dict{String,String}, String}`: A tuple with:
  1. A dictionary mapping each accession (e.g., `"A_1"`) to its species (e.g., `"A"`).
  2. The full path to the output mapping file.

# Examples
```julia-repl
julia> trees = readMultiTopology!("besttrees.tre");

julia> map, path = map_accessions_to_species_dict(trees, "output/", "astral");

julia> map["A_1"]
"A"

julia> println(read(path, String))
A_1 A
B_2 B
```
"""
function map_accessions_to_species_dict(
    gene_trees::Vector{HybridNetwork}, 
    output_dir::String, 
    mode::String = "snaq")

    accession_to_species = Dict{String, String}()

    for tree in gene_trees
        for tip in tree.leaf
            acc = tip.name
            sp = split(acc, "_")[1]
            accession_to_species[acc] = sp
        end
    end

    if mode == "snaq"
        mapping_file = joinpath(output_dir, "id_to_species.csv")
        df = DataFrame(
            individual = String.(keys(accession_to_species)),
            species = String.(values(accession_to_species))
        )
        CSV.write(mapping_file, df)

    elseif mode == "astral"
        mapping_file = joinpath(output_dir, "astral_mapping.txt")
        open(mapping_file, "w") do io
            for (acc, sp) in accession_to_species
                println(io, "$acc $sp")
            end
        end

    else
        error("Invalid mode: $mode. Use \"snaq\" or \"astral\".")
    end

    return accession_to_species, mapping_file
end

"""
    clean_newick_nan(newick_str::String) -> String

Remove NaN (Not-a-Number) values from Newick tree strings.

When ASTRAL outputs trees, it sometimes includes `nan` for missing or undefined
branch lengths (e.g., `A:nan`). These values cause parsing errors in PhyloNetworks.
This function replaces `nan` values with a small default value (0.001) to allow
successful parsing.

# Arguments
- `newick_str::String`: A Newick format tree string potentially containing `nan`

# Returns
- `String`: The cleaned Newick string with `nan` replaced by a smaller number. 

# Examples
```julia-repl
julia> dirty = "((A:nan,B:0.1):0.2,C:0.3);"
julia> clean_newick_nan(dirty)
"((A:0.001,B:0.1):0.2,C:0.3);"
```
"""
function clean_newick_nan(newick_str::String)
    # Replace :nan with :0.000001 (a small default value)
    return replace(newick_str, ":nan" => ":0.000001")
end

"""
    simplify_tip_labels(tree::HybridNetwork) -> HybridNetwork

Simplify tip labels in a hybrid network by removing locus or accession IDs 
from the tip names.

This function modifies each tip label in-place, keeping only the prefix before 
the first underscore. For example, `"A_0"` or `"A_1"` becomes `"A"`.

# Arguments
- `tree::HybridNetwork`: A phylogenetic network whose tip labels will be 
  simplified.

# Returns
- `HybridNetwork`: The same network object, with simplified tip names.

# Examples
```julia-repl
julia> tree = readnewick("(A_0:0.1,B_1:0.2);");

julia> simplify_tip_labels(tree);

julia> tipLabels(tree)
["A", "B"]
```
"""
function simplify_tip_labels(tree::HybridNetwork)
    for tip in tree.leaf
        original = tip.name
        simplified = split(original, "_")[1]  # keep only 'A' from 'A_0'
        tip.name = simplified
    end
    return tree
end


"""
    totallength(net::HybridNetwork, allowmissing::Bool=false)

Sum of all edge lengths in `net`. By default, an error is thrown if an edge
length is missing (coded as -1). Missing edge lengths are skipped if
`allowmissing` is `false`.
Negative (but not -1) edge lengths throw an error log message.
"""
function totallength(net::HybridNetwork, allowmissing::Bool=false)
    tl = 0.0
    for e in net.edge
        el = e.length
        if el == -1.0
            allowmissing && continue # ignore this edge if length is missing
            error("edge $(e.number) does not have a length")
        end
        el < 0 && @error("will add a negative edge length: $(e.length)")
        tl += el
    end
    return tl
end


# The below function is used in filter() to remove rows with repeated species 
# It is retrieved from SNAQ.jl: https://juliaphylo.github.io/SNaQ.jl/dev/man/multiplealleles/ 
"""
    hasrep

Return true if a row (4-taxon set) has a "repeated" species, that is, a species
whose name ends with "__2". Otherwise, return false.

Warning: this function assumes that taxon names are in columns
"t1", "t2", "t3", "t4". For data frames with different column names,
e.g. "taxon1", "taxon2" etc., simply edit the code below by replacing
`:t1` by `:taxon1` (or the appropriate column name in your data).
"""
function hasrep(row)
    occursin(r"__2$", row[:t1]) || occursin(r"__2$", row[:t2]) ||
    occursin(r"__2$", row[:t3]) || occursin(r"__2$", row[:t4])
end 

"""
    set_up_paramname_root(
        dup_rate::Float32, 
        loss_rate::Float32, 
        ratevar::String, 
        n_inds::Int,
        SF::Float32 = 1.0 
    ) -> String 
Generate a standardized parameter name string based based on simulation parameters.
This function creates a concise identifier for a specific set of simulation
parameters, which can be used in file naming and logging.   
"""
function set_up_paramname_root(dup_rate, loss_rate, 
                    ratevar, n_inds, scaling_factor_branch_length, gene_len)
  return "DUP$dup_rate-LOS$loss_rate-RV$ratevar-N_ind$n_inds-SF$scaling_factor_branch_length-genelen$gene_len"
end


"""
    map_tree_based_on_mapping_file(mapping_file::String) 
        -> Dict{String, String} 
Create a mapping of individuals to species from an astral mapping file.
This function reads a mapping file where each line contains an individual
name and its corresponding species name, separated by whitespace. It returns
a dictionary that maps each individual to its species.
"""
function create_mapping_file_based_on_gene_tree(mapping_file::String)
    # Read mapping file and create individual -> species lookup
    mapping_lines = readlines(mapping_file)
    individual_to_species = Dict{String, String}()
    
    for line in mapping_lines
        if !isempty(strip(line))
            parts = split(strip(line))
            if length(parts) == 2
                individual, species = parts
                individual_to_species[individual] = species
            end
        end
    end

    return individual_to_species 
end 

"""
    map_gene_tree_based_on_mapping_file(
        gene_tree::HybridNetwork, 
        individual_to_species::Dict{String, String}
    ) -> HybridNetwork
Map individuals in a gene tree to species using a provided mapping dictionary.
This function takes a gene tree with tips labeled as individuals (e.g., "A_1
or "B_2") and a dictionary mapping individuals to species names (e.g., "A" or "B").
It renames the tips in the gene tree according to the mapping and returns the
modified tree. 
"""
function map_gene_tree_based_on_mapping_file(
    gene_tree::HybridNetwork, 
    individual_to_species::Dict{String, String}) 

    # Map individuals in gene tree to species by renaming tips
    mapped_gene_tree = deepcopy(gene_tree)
    
    # Get tip names and map them
    tip_names = tipLabels(mapped_gene_tree)
    for (i, tip_name) in enumerate(tip_names)
        if haskey(individual_to_species, tip_name)
            species_name = individual_to_species[tip_name]
            # Update tip label
            mapped_gene_tree.leaf[i].name = species_name
        else
            @warn "Individual $tip_name not found in mapping file"
        end
    end

    return mapped_gene_tree 
end 

"""
Prune tips from both trees so that they share the same tip set.
Returns a tuple: (species_tree, gene_tree).
Both trees are deep-copied so the originals remain unchanged.
"""
function prune_nonoverlapping_tips(species_tree::HybridNetwork, 
    gene_tree::HybridNetwork)
    
    species_tree = deepcopy(species_tree)
    gene_tree = deepcopy(gene_tree)
    tiplabels_species_tree = tipLabels(species_tree)
    tiplabels_gene_tree= tipLabels(gene_tree)

    # Tips unique to species_tree
    for tip in setdiff(tiplabels_species_tree, tiplabels_gene_tree)
        deleteleaf!(species_tree, tip)
    end

    # Tips unique to gene_tree
    for tip in setdiff(tiplabels_gene_tree, tiplabels_species_tree)
        deleteleaf!(gene_tree, tip)
    end
    return species_tree, gene_tree
end

"""
    getsisters(node::PhyloNetworks.Node)

Return all sister nodes of the parent edge of `node`.
"""
function get_sister_nodes(node::PhyloNetworks.Node)
    # parent edge & parent node
    e = getparentedge(node)
    p = getparent(node)
    sisters_node = PhyloNetworks.Node[]   
    # loop over all edges attached to parent
    for edge in p.edge
        # skip if it's the same edge
        if edge === e
            continue
        end
        # skip if parent of this edge is not p
        if getparent(edge) !== p
            continue
        end
        push!(sisters_node, getchild(edge))
    end
    return sisters_node
end

"""
    getsisters(leaf::PhyloNetworks.Node) -> Vector{PhyloNetworks.Node}
Return all sister nodes of the parent edge of `leaf`.
Usage: sisterdict = get_sisters_for_all_leafs(tree) 
"""
get_sisters_for_all_leafs(tree::HybridNetwork) =
    Dict(leaf.name => get_sister_nodes(leaf) for leaf in tree.leaf) 

"""
    prune_related_leaves!(
        tree::HybridNetwork, 
        sisterdict::Dict{String, Vector{PhyloNetworks.Node}}
    ) -> pruned_tree::HybridNetwork
"""
function prune_related_leaves!(
    tree::HybridNetwork,
    sisterdict::Dict{String, Vector{PhyloNetworks.Node}}) 

    removed_leaves = Set{String}()

    # iterate over a snapshot of keys 
    for leaf in collect(keys(sisterdict))

        if (leaf in removed_leaves) || !haskey(sisterdict, leaf)
            # if lead_name is pruned before 
            continue
        end 

        leafprefix = split(leaf, "_")[1]
        sisters = sisterdict[leaf] # Vector{Node} 

        for sister in sisters

            if !sister.leaf || isempty(sister.name)
                continue
            end
            sister_prefix = split(sister.name, "_")[1]

            # same prefix prune the sister
            if sister_prefix == leafprefix && !(sister.name in removed_leaves)
                push!(removed_leaves, sister.name)
                # remove from the network (ok if it's already gone)
                try
                    deleteleaf!(tree, sister.name)
                catch
                    @warn "Failed to delete leaf $sister.name from tree"
                end

                # drop as a dict key if present
                if haskey(sisterdict, sister.name)
                    delete!(sisterdict, sister.name)
                end
            end
        end
    end
    return tree 
end 


"""
Iterates over a set of examples.

# Description
This loop is intended to process or analyze each example in a collection. 
Replace `example` with the specific variable or data structure representing your examples.

# Example
for example: 
    dict = group_tips_by_species(net) 
Dict{String, Vector{String}} with 8 entries:
"B" => ["B_0"]
"A" => ["A_1"]
"C" => ["C_1", "C_0"]
"D" => ["D_1"]
"G" => ["G_1", "G_0"]
"E" => ["E_0", "E_1"]
"H" => ["H_1", "H_0"]
"F" => ["F_1"]
"""
function group_tips_by_species(tree)
        species_dict = Dict{String, Vector{String}}()
        for tip in tipLabels(tree)
            sp = split(tip, "_")[1]
            push!(get!(species_dict, sp, String[]), tip)
        end
        return species_dict
    end  

"""
    all_combinations(species_dict::Dict{String, Vector{String}}) 
        -> Iterators.ProductIterator 
Generate all combinations of tip selections from a species-to-tips dictionary.
This function takes a dictionary where keys are species names and values are
vectors of tip labels belonging to those species. It returns an iterator that
produces all possible combinations of selecting one tip from each species.
Example: 
The below example will have 1 x 1 x 2 x 1 x 2 x 2 x 2 x 1 = 16 combinations 
```julia-repl 
dict = group_tips_by_species(net)  
Dict{String, Vector{String}} with 8 entries:
  "B" => ["B_0"]
  "A" => ["A_1"]
  "C" => ["C_1", "C_0"]
  "D" => ["D_1"]
  "G" => ["G_1", "G_0"]
  "E" => ["E_0", "E_1"]
  "H" => ["H_1", "H_0"]
  "F" => ["F_1"]
julia> for combo in all_combinations(dict)
           println(collect(combo))
       end 
# omit the printing here
length(collect(all_combinations(dict))) 
16 
```
"""
function all_combinations(species_dict)
    choices = [v for v in values(species_dict)]
    return Iterators.product(choices...)
end 

"""
    subtree_for_combination(
        tree::HybridNetwork, 
        tips_to_keep::Vector{String}
    ) -> HybridNetwork 
Create a subtree by keeping only specified tips from the original tree.
This function takes a phylogenetic tree and a list of tip labels to retain.
It returns a new tree that includes only those tips, pruning all others.    
"""
function subtree_for_combination(tree, tips_to_keep)
    copy_tree = deepcopy(tree)
    tips_to_drop = setdiff(tipLabels(copy_tree), tips_to_keep)
    for t in tips_to_drop
        deleteleaf!(copy_tree, t) # or your preferred prune function
    end
    return copy_tree
end

"""
    scale_branch_lengths(
        net::HybridNetwork, 
        scale_factor::Float64
    ) -> HybridNetwork

Return a new hybrid network with all branch lengths scaled by a given factor.

This function creates a deep copy of the input network and multiplies the length
of each edge by the specified scale factor. The original network is not modified.

# Arguments
- `net::HybridNetwork`: The hybrid network whose branch lengths will be scaled.
- `scale_factor::Float64`: The factor by which to scale each branch length.

# Returns
- `HybridNetwork`: A new network with scaled branch lengths.

# Examples
"""
function scale_branch_lengths(
    net::HybridNetwork, 
    scale_factor::Float64
    ) 
    modified_net = deepcopy(net)
    for edge in modified_net.edge
        edge.length *= scale_factor 
    end 
    return modified_net 
end

"""
    scale_species_tree_with_sub_scaling(
        species_tree_with_sub_scaling::String, 
        scaling_factor::Float64
    ) -> String

Scale all branch lengths (numbers before "*") in a species tree string 
with substitution scaling by a given factor.

This function processes a species tree string where branch lengths are 
formatted as "length*rate" and multiplies each length component by the 
specified scaling factor, preserving the rate component.

# Arguments
- `species_tree_with_sub_scaling::String`: A Newick tree string with branch 
  lengths in the format "length*rate" (e.g., "A:3440*1.778")
- `scaling_factor::Float64`: The factor by which to scale each branch length

# Returns
- `String`: A modified Newick string with scaled branch lengths, preserving 
  the rate components

# Examples
```julia-repl
julia> tree = "(A:3440*1.778,B:880*0.189);"
julia> scale_species_tree_with_sub_scaling(tree, 2.0)
"(A:6880.0*1.778,B:1760.0*0.189);"
```
"""
function scale_species_tree_with_sub_scaling(
    species_tree_with_sub_scaling::String, 
    scaling_factor::Float64
    ) 
    # Pattern to match branch lengths in the format ":number*"
    pattern = r":(\d*\.?\d+(?:[eE][+-]?\d+)?)\*" 
    
    # Replace each match by scaling the number before "*"
    result = replace(species_tree_with_sub_scaling, pattern => function(match)
        # Extract the number part (without ":" and "*")
        number_str = match[2:end-1]  # Remove ":" at start and "*" at end
        number = parse(Float64, number_str)
        scaled_number = number * scaling_factor
        return ":$(scaled_number)*"
    end)
    
    return result
end 

"""
    modify_locus_tree_labels!(locus_tree::HybridNetwork) -> HybridNetwork
Modify tip labels in a locus tree by removing anything after "_"
"""
function modify_locus_tree_labels!(locus_tree::HybridNetwork) 
    for tip in locus_tree.leaf
        original = tip.name
        modified = split(original, "_")[1]  # take the part before the first underscore
        tip.name = modified
    end
    return locus_tree
end 

"""
    get_alignment_length(fasta_file::String) -> Int

Get the sequence length of the first sequence in a FASTA alignment file.

Since all sequences in a proper alignment have the same length, this function
reads only the first sequence and returns its length.

# Arguments
- `fasta_file::String`: Path to a FASTA format alignment file.

# Returns
- `Int`: The sequence length in base pairs (including gaps).

# Examples
```julia-repl
julia> get_alignment_length("alignment.fasta")
30  # First sequence has 30 bp
```

# Notes
- Handles multi-line FASTA sequences (sequence split across multiple lines)
- Strips whitespace from sequence lines
- Gap characters ('-') are counted as part of the sequence length
- Stops reading after the first complete sequence for efficiency
"""
function get_alignment_length(fasta_file::String)
    open(fasta_file, "r") do f
        current_seq = ""
        started = false
        
        for line in eachline(f)
            if startswith(line, ">")
                if started && !isempty(current_seq)
                    # Found the next header, so we're done with first sequence
                    return length(current_seq)
                end
                started = true
            elseif started
                # Accumulate sequence (removing whitespace)
                current_seq *= strip(line)
            end
        end
        
        # Return length of the first sequence
        return length(current_seq)
    end
end 

"""
    diagnose_simphy_locus_trees(simphy_output_dir::String, rep::Int, iter::Int) -> String

Diagnose issues with locus trees generated by SimPhy by examining the l_trees.trees file.

This function checks for common issues that cause SimPhy failures:
1. Whether the l_trees.trees file was generated
2. Whether any locus tree has only taxa from a single species after removing "Lost-" taxa
3. Whether any locus tree has insufficient species diversity

# Arguments
- `simphy_output_dir::String`: Path to the SimPhy output directory for this iteration
- `rep::Int`: Replicate number (for error messages)
- `iter::Int`: Iteration number (for error messages)

# Returns
- `String`: Detailed diagnostic message describing the failure

# Examples
```julia-repl
julia> diagnose_simphy_locus_trees("/path/to/Int5/1", 1, 5)
"SimPhy generated 50 locus trees but 2 are problematic:
  Locus tree #3: Only 1 taxon (A) remains after removing Lost taxa
  Locus tree #17: All 2 remaining taxa are from species 'H'"
```

# Notes
- Searches for l_trees.trees in both output_dir and output_dir/1
- Removes all tips starting with "Lost-" before analysis
- Extracts species prefix by taking the part before the first underscore
- Reports specific locus tree line numbers (1-indexed)
"""
function diagnose_simphy_locus_trees(simphy_output_dir::String, rep::Int, iter::Int)
    # Check if l_trees.trees file exists
    # SimPhy may create output in a subdirectory named "1"
    l_trees_path = joinpath(simphy_output_dir, "l_trees.trees")
    if !isfile(l_trees_path)
        # Try the "1" subdirectory
        l_trees_path_alt = joinpath(simphy_output_dir, "1", "l_trees.trees")
        if !isfile(l_trees_path_alt)
            return "SimPhy failed before generating locus tree file (l_trees.trees not found)"
        end
        l_trees_path = l_trees_path_alt
    end
    
    # Read all locus trees
    locus_trees = readlines(l_trees_path)
    num_locus_trees = length(locus_trees)
    
    # Track problematic trees
    problematic_trees = []
    
    # Analyze each locus tree
    for (locus_id, newick_str) in enumerate(locus_trees)
        if isempty(strip(newick_str))
            push!(problematic_trees, (locus_id, "Empty tree string"))
            continue
        end
        
        try
            # Parse the Newick string into a HybridNetwork
            net = readnewick(newick_str)
            
            # Get all tip names
            all_tips = [tip.name for tip in net.leaf]
            
            # Identify tips to remove (starting with "Lost-")
            lost_tips = filter(name -> startswith(name, "Lost-"), all_tips)
            
            # Remove "Lost-" tips
            for lost_tip in lost_tips
                try
                    deleteleaf!(net, lost_tip)
                catch e
                    # If deletion fails, note it but continue
                    push!(problematic_trees, (locus_id, "Failed to remove Lost tip: $lost_tip"))
                    continue
                end
            end
            
            # Get remaining tips after removing "Lost-" taxa
            remaining_tips = [tip.name for tip in net.leaf]
            num_remaining = length(remaining_tips)
            
            # Extract species prefixes (part before first underscore or the whole name)
            species_list = []
            for tip_name in remaining_tips
                # Split by underscore and take first part as species
                parts = split(tip_name, "_")
                species = parts[1]
                push!(species_list, species)
            end
            
            # Get unique species
            unique_species = unique(species_list)
            num_unique_species = length(unique_species)
            
            # Check for problematic patterns
            if num_remaining == 0
                push!(problematic_trees, (locus_id, "All taxa were Lost (no remaining taxa)"))
            elseif num_remaining == 1
                # Only one taxon remains
                species = unique_species[1]
                push!(problematic_trees, (locus_id, "Only 1 taxon ($species) remains after removing Lost taxa"))
            elseif num_unique_species == 1
                # Multiple taxa but all from same species
                species = unique_species[1]
                push!(problematic_trees, (locus_id, "All $num_remaining remaining taxa are from species '$species' ($(join(remaining_tips, ", ")))"))
            elseif num_unique_species < 4
                # Less than 4 species (not phylogenetically informative for 8-taxon tree)
                push!(problematic_trees, (locus_id, "Only $num_unique_species species remain: $(join(unique_species, ", ")) (need ≥4 for phylogenetic inference)"))
            end
            
        catch e
            # Failed to parse or process the tree
            push!(problematic_trees, (locus_id, "Failed to parse or process tree: $e"))
        end
    end
    
    # Generate diagnostic message
    if isempty(problematic_trees)
        return "SimPhy failed after generating $num_locus_trees locus trees, but cause is unknown. All locus trees appear valid. Inspect the locus tree file: $l_trees_path"
    else
        msg = "SimPhy generated $num_locus_trees locus trees but $(length(problematic_trees)) are problematic:\n"
        for (locus_id, reason) in problematic_trees
            msg *= "  Locus tree #$locus_id: $reason\n"
        end
        msg *= "Locus tree file: $l_trees_path"
        return msg
    end
end