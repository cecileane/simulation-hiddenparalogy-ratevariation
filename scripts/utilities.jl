#= Utility functions shared across all pipeline scripts.

Part 0 -- Logging helpers for distributed workers
    worker_println  -> println that also writes to worker log file
    open_worker_log -> open per-worker log at start of parallel job
    close_worker_log -> close per-worker log at end of parallel job

Part 1 -- Gene tree cleanup after SimPhy
    replace_taxa_name   -> rename one tip given species char + newick
    count_missing       -> count missing/duplicate taxa in a newick
    modify_newicks       -> filter and simplify a gene tree file
    modify_newick_for_n_genes -> run modify_newicks across n gene files

Part 2 -- Simulation control helpers
    seed_generator  -> build N x M seed arrays from a master seed
    calculate_batch  -> estimate how many SimPhy iters remain

Part 3 -- General utilities
    pad_number           -> zero-padded string matching SimPhy output
    check_existing_dir   -> interactive prompt to clear output dirs
    replace_tips_with_letters -> convert species names to letters
    replace_pop_inInd_file -> assign unique pop IDs in .ind files
    generate_software_seeds -> per-software seeds from master seed
    generate_master_seed   -> deterministic seed from param dict
    get_dict_for_seed_setting -> extract seed params from name string
    totallength     -> sum all branch lengths in a tree/network
    hasrep          -> detect repeated taxa in a newick string
    set_up_paramname_root  -> build the canonical param name string
=#
using PhyloNetworks
using StatsBase
using StableRNGs
using Primes
using Statistics
using IterTools
using Glob
using Distributions

#-----------------------------------------------#       
#               Part 0  
#    Logging helper functions
#-----------------------------------------------#
"""
println that also writes to the log file of this worker, if logging is on
Input: anything println takes
Notes: uses the globals log_output (Bool) and worker_logfile (IOStream or nothing)
"""
function worker_println(args...)
    println(args...)
    flush(stdout)
    log_active = @isdefined(log_output) && log_output
    file_open  = @isdefined(worker_logfile) && worker_logfile !== nothing
    if log_active && file_open
        try
            println(worker_logfile, args...)
            flush(worker_logfile)
        catch e
            # ignore, a logging problem should not kill the run
        end
    end
end

"""
Open the log file of this worker: screen_<paramname>_worker<id>.log in outfolder
Call it at the start of each worker. 
Uses the globals log_output, outfolder, paramname_root
"""
function open_worker_log()
    if @isdefined(log_output) && log_output
        try
            log_filename = joinpath(outfolder,
                "screen_$(paramname_root)_worker$(myid()).log")
            global worker_logfile = open(log_filename, "w")
            println(worker_logfile, "=== Worker $(myid()) Log Started ===")
            flush(worker_logfile)
        catch e
            @warn "Could not open worker log file: $e"
        end
    end
end

"""
Close the log file opened by open_worker_log, at the end of each worker
"""
function close_worker_log()
    if @isdefined(worker_logfile) && worker_logfile !== nothing
        try
            println(worker_logfile, "=== Worker $(myid()) Log Ended ===")
            close(worker_logfile)
            global worker_logfile = nothing
        catch e
            # ignore
        end
    end
end

#-----------------------------------------------#       
#               Part 1  
#    Modify newick strings generated from Simphy
#    Goal: Mimic hidden paralogy    
#-----------------------------------------------#
"""
Inputs:
    -taxon (Char): 
    The target character to search for in tip labels.
    -newick (String): 
    The Newick tree string containing tip labels in the format char_x_y.
Outputs: Tuple (count, modified_tree):
    -count (Integer): Number of occurrences of char in the tree.
    -modified_tree (String): The Newick string, modified if applicable.
    1) If tip appears in the string for only once 
        --> Return (count, modified tree).
    2) If tip appears in the string for 0 or more than once 
        --> Return (count, unmodified tree).
Examples:
# Case 1: Char appears in the string for only once
Input: replace_taxa_name('A', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
Output: (1, "(A_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")

# Case 2: Char appears in the string for 0 time
Input: replace_taxa_name('F', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
Output: (0, "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")

# Case 3: Char appears in the string for >= 2 times:
Input: replace_taxa_name('A', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_1_0:0.02,A_2_1:0.02):0.01)));")
Output: (2, "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_1_0:0.02,A_2_1:0.02):0.01)));")
# Note there are A_0_1, A_1_0, A_2_1. It outputs 2 here because A_0_1 and A_2_1 are gene copies from the same accession.

# Case 4: Similar to above, char appears >= 2 times, but these are different individuals from the same species:
Input: replace_taxa_name('A', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_0_0:0.02,A_0_2:0.02):0.01)));")
Output: (1, "(A_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_0:0.02,A_2:0.02):0.01)));")
because all three As are gene copies from different accessions.
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
    # max copies of this taxon across accessions;
    # 0 = absent, >=2 = duplicated gene (non-hidden paralogy)
    if num == 0 || num >= 2
        return (num, newick)
    else # If # of taxons == 1, remove locus ID from the tip 
        modified_newick = replace(newick, regex => s"\1_\2") 
        return(num, modified_newick)
    end
end

"""
Inputs:
    -newick (String): A Newick tree string containing taxa in the format char_x_y. 
    Taxa are hardcoded as "ABCDEFGH".
Outputs:
    -Tuple (missing_count, modified_tree) if all gene copies are unique:
        missing_count (Integer): Number of taxa missing from the tree.
        modified_tree (String): Modified Newick string with simplified taxa labels.
    -Returns nothing if there are repeated gene copies 
        (then the paralogy is not hidden).

Examples:
# Case 1: All taxa appear once but F, G, H are missing (number of missing = 3).
Input: count_missing("(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
Output: (3, "(A_1:0.01,(((B_1:0.004,C_1:0.003):0.007,(D_0:0.02,E_1:0.02):0.01)));")

# Case 2: B, G, H missing (number of missing = 3) and there are repeated copies A_0_1 and A_1_1. Returns nothing
Input: count_missing("(F_0_1:0.8,(A_0_1:0.01,((A_1_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
Output: nothing
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
    return (num_missing, nwk)
end

"""
Input: one file with one newick tree from simphy, and the names of the output file and of the
    intermediate .mapsl / .maplg files of that gene 
    (removed with the tree if the tree is dropped)
Output: the modified newick is written to the output file, 
    returns (n_repeated, n_missing)
- If the input newick has repeated taxa,
    not included in the downstream analysis (n_repeated = 1)
- If the input newick has more than max_taxa_missing taxa missing, 
    not included (n_missing = 1)
Notes: If the input has repeated taxa, then paralogy is not hidden. 
    To simulate hidden paralogy, the tree is not processed.
    If no hidden paralogy got simulated, 
    then dup_rate = 0 and loss_rate = 0, so no gene tree is removed.
However, this code still renames the tree tips so that A_locusID_accessionID 
is changed to A_accessionID.
"""
function modify_newicks(input::String, output::String, 
                        intermediate_file_mapsl::String, 
                        intermediate_file_maplg::String,
                        max_taxa_missing::Int,
                        debug_mode::Bool=false)
    if !isfile(input)
        @warn "Input file $input does not exist, skipping"
        return (0, 0)
    end
    
    # read the whole file at once: on AFS, open/read/close can interleave with a concurrent write
    # one simphy file = one newick
    nwk = try
        content = read(input, String)
        String(strip(content))
    catch e
        @error "Failed to read $input" exception=(e, catch_backtrace())
        return (0, 0)
    end
    
    if isempty(nwk)
        @warn "Input file $input is empty, skipping"
        return (0, 0)
    end
    
    tree_repeated_taxa = 0 
    tree_insufficient_taxa = 0 
    
    temp = try
        count_missing(nwk)
    catch e
        @error "Failed to process tree in $input" exception=(e,
            catch_backtrace())
        return (0, 0)
    end

    # remove a file if it exists, warn if that fails
    function try_rm(path)
        isempty(path) && return
        try
            isfile(path) && rm(path; force=true)
        catch e
            @warn "Failed to remove $path" exception=(e, catch_backtrace())
        end
    end

    if isnothing(temp) # repeated taxa
        println("The input $input has repeated taxa")

        if !debug_mode
            try_rm(intermediate_file_mapsl)
            try_rm(intermediate_file_maplg)
            try_rm(input)
        else
            println("DEBUG: Keeping intermediate files and tree for analysis")
        end

        tree_repeated_taxa = 1

    elseif temp[1] > max_taxa_missing
        println("The input $input has too many taxa missing " *
            "(missing=$(temp[1]), max=$max_taxa_missing)")

        if !debug_mode
            try_rm(intermediate_file_mapsl)
            try_rm(intermediate_file_maplg)
            try_rm(input)
        else
            println("DEBUG: Keeping intermediate files and tree for analysis")
            println("DEBUG: Tree content: $nwk")
            println("DEBUG: Modified would be: $(temp[2])")
        end 

        tree_insufficient_taxa = 1
        
    else # good tree, write it
        modified_tree = temp[2]
        try
            open(output, "w") do io
                write(io, modified_tree * "\n")
            end
        catch e
            @warn "Failed to write $output" exception=(e, catch_backtrace())
        end
    end
    
    return tree_repeated_taxa, tree_insufficient_taxa  
end


"""
Run modify_newicks on the n gene trees (newicks) of one simphy iteration
Input:
    -input: input dir path, with g_trees<geneID>.trees (zero-padded IDs)
    -output: output dir path, gets g_trees_noLocusID_Gene<geneID>_Int<iteration>.trees
    -n: Number of gene trees / newicks to be modified
    -iteration, dup_rate, n_inds, max_taxa_missing, debug_mode: see simulation.jl
Output: see the return block at the end of the function (counts of dropped and kept trees by category)
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
    
    # counts of filtered trees (removed from dataset)
    tree_repeated_taxa_tot = 0
    tree_insufficient_taxa_tot = 0

    # counts of kept trees by event type
    num_trees_experiencing_gene_duplication_and_loss = 0
    num_trees_experiencing_gene_loss_only = 0
    # trees with no duplication and no loss
    num_trees_experiencing_nothing = 0

    # leaf counts across kept trees (for average reporting)
    num_leaf_left_in_trees = Int[]
    
    # Save a list of gene trees name for each category 
    gene_trees_gene_duplication_and_loss = String[]
    gene_trees_gene_loss_only = String[]
    # No need to save for gene tree experiencing nothing
    # because the rest are gene trees experiencing nothing

    for gene_tree in 1:n
        genenum_string = pad_number(gene_tree, n)
        output_tree_name = "g_trees_noLocusID_Gene$(genenum_string)" *
            "_Int$(iteration).trees"
        input_dir = joinpath(input, "g_trees$genenum_string.trees")
        output_dir = joinpath(output, output_tree_name)

        # Keep intermediate .mapsl/.maplg files only for valid trees
        # so postprocessing can classify loss-only vs dup+loss events.
        intermediate_file_mapsl = joinpath(input, "$genenum_string.mapsl")
        intermediate_file_maplg = joinpath(input, "$(genenum_string)l1g.maplg")

        tree_repeated_taxa, tree_insufficient_taxa = modify_newicks(
            input_dir, output_dir,
            intermediate_file_mapsl, intermediate_file_maplg,
            max_taxa_missing, debug_mode)

        tree_repeated_taxa_tot += tree_repeated_taxa
        tree_insufficient_taxa_tot += tree_insufficient_taxa

        # Classify kept trees by event type using the .mapsl file:
        #   loss only  -> no "dup" in file, num leaf < 8*n_inds
        #   dup + loss -> "dup" present in file
        #   nothing    -> no dup, no loss (all 8*n_inds leaves kept)
        # Hidden-paralogy classification is deferred to postprocessing.
        
        # If dup_rate == 0, no mapsl file so skip the below analysis 
        if dup_rate > 0 
            # When duplicate >= 1 and/or loss >= 1: 
            # Both mapsl and maplg files should exist 
            # When duplicate == 0 and loss >= 0: 
            # only maplg file should exist 
            if !isempty(intermediate_file_mapsl) &&
                    isfile(intermediate_file_mapsl)
                lines = lowercase.(readlines(intermediate_file_mapsl))
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
    avg_num_leaf_left_in_trees = isempty(num_leaf_left_in_trees) ?
        NaN : mean(num_leaf_left_in_trees)
    println("Avg leaves per kept tree: $avg_num_leaf_left_in_trees")

    #------ Finally, clean up intermediate files -------# 
    # remove intermediate files to save space under non-debugging mode  
    if debug_mode == false # debugging mode will keep intermediate files 
        
        # search maplg only: when dup=0, only maplg exists (no mapsl)
        if !isempty(collect(glob("*.maplg", input)))
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
            # no mapsl or maplg file: remove the iteration folder
            # one level up from the input directory
            parent_dir = dirname(input) # This is the interation folder 
            if isdir(parent_dir) 
                rm(parent_dir; force=true, recursive=true)
            end 
        end 

        # move the .trees files one level up, to the iteration folder, and remove the nested Int/1 folder
        if isdir(input)
            iteration_dir = dirname(input)
            trees_files = filter(f -> endswith(f, ".trees"), readdir(input))
            for file in trees_files
                src_path = joinpath(input, file)
                dst_path = joinpath(iteration_dir, file)
                mv(src_path, dst_path; force=true)
            end
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

"""
    blank_unused_locus_trees(simphyfolder, singlecopyfolder)

In each Int<i>/l_trees.trees (or Int<i>/1/l_trees.trees) of one replicate,
keep only the locus trees of the kept gene trees and blank all other lines.
A kept gene g_trees_noLocusID_Gene<g>_Int<i>.trees uses line <g>, which is
what simulation_postprocess.jl reads, so line numbers must not change.
The locus trees of discarded gene families are never used and take most of
the disk space when dup/loss rates are high.
"""
function blank_unused_locus_trees(simphyfolder::String,
                                  singlecopyfolder::String)
    # kept gene ids for each iteration
    kept = Dict{Int, Set{Int}}()
    for file in readdir(singlecopyfolder)
        m = match(r"^g_trees_noLocusID_Gene(\d+)_Int(\d+)\.trees$", file)
        m === nothing && continue
        push!(get!(kept, parse(Int, m.captures[2]), Set{Int}()),
              parse(Int, m.captures[1]))
    end

    for int_folder in readdir(simphyfolder)
        m = match(r"^Int(\d+)$", int_folder)
        m === nothing && continue
        # legacy structure has a "1" subfolder, as in simulation_postprocess.jl
        temp_dir = joinpath(simphyfolder, int_folder, "1")
        locus_tree_file = isdir(temp_dir) ?
            joinpath(temp_dir, "l_trees.trees") :
            joinpath(simphyfolder, int_folder, "l_trees.trees")
        isfile(locus_tree_file) || continue

        genes = get(kept, parse(Int, m.captures[1]), Set{Int}())
        locus_trees = readlines(locus_tree_file)
        for g in genes # every kept gene must have its locus tree
            if g > length(locus_trees) || isempty(locus_trees[g])
                error("No locus tree for gene $g in $locus_tree_file")
            end
        end
        new_locus_trees = [i in genes ? t : ""
                           for (i, t) in enumerate(locus_trees)]
        new_locus_trees == locus_trees && continue # nothing to blank

        temp_file = locus_tree_file * ".tmp"
        open(temp_file, "w") do io
            for t in new_locus_trees
                println(io, t)
            end
        end
        mv(temp_file, locus_tree_file; force=true)
    end
end

#-----------------------------------------------#
#               Part 2
#    Re-run simphy until:
# 1) n_genes hit the target; or 
# 2) max_iteration reaches   
#-----------------------------------------------#
"""
Writes n x m random integers to be used as seeds
Input:  -master_seed: An Int for StableRNG
        -n: number of rows (reps)
        -m: number of columns (seeds per rep)
        -output_dir, output_file_name: where to write the seeds 
        (tab separated, header col1 ... colm)
Output: - the seed file in output_dir
        - Returns the n x m array of seeds (Int32, positive)
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
Estimates the batch size for the next simulation based on missing rate and target gene count.
Inputs:
    - num (int): Number of gene trees existing
    - n_genes (int): Total number of target gene trees.
Output:
    -int: Rounded estimated batch size for the next simulation.
Calculation explanations:
    For a missing rate of 0.2 and 30 target genes:
    - Missing trees: 30 * 0.2 = 6
    - Success rate: 1 - 0.2 = 0.8
    - estimated batch size to be simulated: 6 / 0.8 = 7.5 -> 8
    If no tree succeeded yet, return n_genes. If the success rate is <= 0.2, cap at 6 * n_genes
"""
function calculate_batch(num::Int, n_genes::Int) 
    num_missing = n_genes - num 
    success_rate = num / n_genes 
    if success_rate == 0
        estimated_batch = n_genes 
    elseif success_rate <= 0.2
        println("success rate too small, estimated batch size = 5 * n_genes")
        estimated_batch = n_genes * 6
        # cap batch size; if success_rate is tiny the raw estimate explodes
    else 
        estimated_batch = round(num_missing / success_rate) 
    end 
    # println("For replicate $n_reps estimated batch size = $estimated_batch")
    return Int(estimated_batch)
end 


#-----------------------------------------------#       
#               Part 3  
#       Other utilities functions     
#-----------------------------------------------#
"""
Converts num into a zero-padded string based on the number of digits required
by range. This function is required because simphy outputs numbers as below.
Inputs:
    - num::Int: The number to pad.
    - range::Int: The range determining the padding width.
Outputs:
    - String: A zero-padded string representation of num.
Examples:
julia> pad_number(1, 10)   # "01"
julia> pad_number(1, 5)    # "1"
julia> pad_number(1, 200)  # "001"
"""
function pad_number(num::Int, range::Int) 
  digits = ceil(Int, log10(range + 1))
  number_string = lpad(string(num), digits, '0') 
  return number_string
end

"""
Which taxa are the hybrid, and who are the donors, in an H=1 network
Input:
    -n (HybridNetwork): network with 1 hybrid node. It gets rerooted in place, 
    pass a deepcopy to keep the original
    -outgroup (String): taxon to root at, default "A"
Output: NamedTuple (hybrid_taxon, major_donor, minor_donor), 
    each a comma separated list of taxa, e.g. "B,C"
    "NA" for all three if something goes wrong (not exactly 1 hybrid, rooting failed, ...)
Notes: if rootatnode! fails with RootMismatch (outgroup is below the hybrid), 
    root on the major hybrid edge instead
Example: info = get_hybrid_info(net, "A"); info.hybrid_taxon  # "G"
"""
function get_hybrid_info(n::HybridNetwork, outgroup::String="A")

    # leaves below start_node, going down only and never through stop_node
    function leaves_below(start_node, stop_node=nothing)
        leaves = String[]
        stack = [start_node]
        visited = Set{Int}()
        while !isempty(stack)
            curr = pop!(stack)
            curr.number in visited && continue
            push!(visited, curr.number)
            if curr.leaf
                push!(leaves, curr.name)
            else
                for e in curr.edge
                    c = getchild(e)
                    if getparent(e) === curr && c !== stop_node
                        push!(stack, c)
                    end
                end
            end
        end
        return sort(leaves)
    end

    # donor clade of one hybrid edge: leaves below its parent node, 
    # not going through the hybrid node
    # if that is empty (the parent's only child is the hybrid node), 
    # climb up along tree edges until some leaves are reachable. 
    # An undirected traversal would spill into the whole network.
    # e.g. in (((F:0.144,((D:0.077,E:0.075)Rrrl:0.045,((B:0.04,C:0.04)B_a:0.054,
    #   (#H1:0.0::0.054)admixa:0.0)G_z:0.004)Rrrz:0.012)Rrr:0.021,
    #   (H:0.159,((G:0.081)#H1:0.0::0.946)admixb:0.072)rooty:0.134)root:0.126,A:0.126)rootb;
    #   admixb has #H1 as only child -> climb to rooty, leaves_below(rooty, #H1) = {H}, the major donor
    #   admixa has #H1 as only child -> climb to G_z, leaves_below(G_z, #H1) = {B, C}, the minor donor
    function donor_leaves(start_node, hnode)
        current = start_node
        while true
            lv = leaves_below(current, hnode)
            isempty(lv) || return lv

            # climb along the tree edge only, not along a hybrid edge
            tree_parent_edge = nothing
            for e in current.edge
                if getchild(e) === current && !e.hybrid
                    tree_parent_edge = e
                    break
                end
            end
            isnothing(tree_parent_edge) && return String[]  # hit root
            current = getparent(tree_parent_edge)
        end
    end

    try
        # root on the outgroup
        try
            rootatnode!(n, outgroup)
        catch e
            if isa(e, PhyloNetworks.RootMismatch)
                target_hnode = nothing
                for hnode in n.hybrid
                    if outgroup in leaves_below(hnode)
                        target_hnode = hnode
                        break
                    end
                end
                if isnothing(target_hnode)
                    error("get_hybrid_info: could not find hybrid " *
                          "containing '$outgroup' as descendant.")
                end
                e_maj, _, _ = PhyloNetworks.hybridEdges(target_hnode)
                rootonedge!(n, e_maj)
            else
                rethrow(e)
            end
        end

        directedges!(n)

        if length(n.hybrid) != 1
            @warn "get_hybrid_info: expected 1 hybrid node, " *
                  "found $(length(n.hybrid)). Returning NA."
            return (hybrid_taxon="NA", major_donor="NA", minor_donor="NA")
        end

        hnode   = n.hybrid[1]
        e_major, e_minor, _ = PhyloNetworks.hybridEdges(hnode)

        hybrid_leaves = leaves_below(hnode)
        major_leaves  = donor_leaves(getparent(e_major), hnode)
        minor_leaves  = donor_leaves(getparent(e_minor), hnode)

        join_leaves(v) = isempty(v) ? "NA" : join(v, ",")

        return (hybrid_taxon = join_leaves(hybrid_leaves),
                major_donor  = join_leaves(major_leaves),
                minor_donor  = join_leaves(minor_leaves))

    catch e
        @warn "get_hybrid_info: unexpected error, returning NA. Error: $e"
        return (hybrid_taxon="NA", major_donor="NA", minor_donor="NA")
    end
end

"""
check_existing_dir(dir_input::Vector{String})

Prompt the user to delete existing output directories before a run.
Useful when re-running simulations and you want a clean slate.

Prompts: ALL (remove all at once), Y/N (per directory).
Throws an error if the user declines removal of an existing directory.
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
                println("Ok, removing each directory individually...")
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
                println("Directory $dir exists. Remove it? Y/N")
                user_input = lowercase(readline()) 
                if user_input == "y"
                    rm(dir; recursive=true)
                    break
                elseif user_input == "n"
                    error("$dir already exists and removal was declined")
                    return
                else
                    println("Invalid input. Please type Y or N.")
                end
            end
        end
    end
end

"""
This function replaces text tips with letters (max # taxa: 26). Examples:

1. String input -> string output:
Input: replace_tips_with_letters("(Homo:0.01,(((Cat:0.004,Dog:0.003):0.007,(Dinosaur:0.02,Birds:0.02):0.01)));")
Output: "(A:0.01,(((B:0.004,C:0.003):0.007,(D:0.02,E:0.02):0.01)));"

2. Tree input -> tree output:
Input: replace_tips_with_letters(readnewick("(Homo:0.01,(((Cat:0.004,Dog:0.003):0.007,(Dinosaur:0.02,Birds:0.02):0.01)));"))
Output:
HybridNetwork, Rooted Network
9 edges, 10 nodes (5 tips, 5 internal).
Tip labels: A, B, C, D, ...
"(A:0.01,(((B:0.004,C:0.003):0.007,(D:0.02,E:0.02):0.01)));"

3. Handles special characters like *:
Input: replace_tips_with_letters("(Homo:0.01*4,(((Cat:0.004*5,Dog:0.003*10):0.007,(Dinosaur:0.02*20,Birds:0.02*3):0.01*9)));")
Output: "(A:0.01*4,(((B:0.004*5,C:0.003*10):0.007,(D:0.02*20,E:0.02*3):0.01*9)));"
If there are more than 26 tips, returns the input with a warning.
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
            println("Warning: more than 26 tips, skipping label replacement.")
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
            println("Warning: more than 26 tips, skipping label replacement.")
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
Input:
    -folder_path_list (Vector): folder of each replicate
    -simulation_rep (Int): which replicate
    -path_string (String): subfolder or file name to add
Output: joinpath(folder_path_list[simulation_rep], path_string)
e.g. setup_rep_output_folders(["rep1_dir", "rep2_dir"], 2, "snaqfolder") gives "rep2_dir/snaqfolder"
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
    replace_pop_inInd_file(input_file::String, output_file::String)

Replace the "POP" placeholder in column 3 of an eigenstrat `.ind` file
with the species letter extracted from column 1 (e.g. "A_0" -> "A").
Input must be tab-separated with exactly three columns per line.
"""
function replace_pop_inInd_file(input_file::String, output_file::String)
    open(input_file, "r") do infile
        open(output_file, "w") do outfile
            for line in eachline(infile)
                fields = split(line, '\t')
                if length(fields) != 3
                    error("$(input_file): .ind must have three columns")
                end 
                if fields[3] != "POP"
                    error("$(input_file): column 3 must be 'POP'")
                end 
                fields[3] = split(fields[1], '_')[1]
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

Map each software name to a stable integer seed derived from `master_seed`.
Uses StableRNG so the same inputs always produce the same seeds across
Julia versions.

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
function generate_software_seeds(master_seed::Int,
        software_names::Vector{String})
    rng = StableRNG(master_seed)
    seeds = abs.(rand(rng, Int, length(software_names)))
    return Dict(
        software_names[i] => seeds[i] for i in eachindex(software_names))
end

"""
    generate_master_seed(params::Dict{String, Any}) -> Int

Build a deterministic integer seed from simulation parameters by encoding
values as products of primes. Floats use per-digit primes; categoricals
(e.g. ratevar) use a fixed lookup; integers map to the next prime.
Uses `big()` internally to avoid overflow.
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
            # nextprime: smallest prime >= value
        else
            error("""Oh no something is wrong with the parameters setting. 
                    Unable to generate master seed""")
        end   
        master_seed *= prime
    end
    return Int(master_seed) 
    # generate_software_seeds requires Int
end

"""
    sample_substitution_params(seed)

Draw per-gene HKY substitution parameters from stable distributions:
- kappa    ~ LogNormal(1.4215, 0.2798)
- basefreqs ~ Dirichlet(66.59, 38.41, 38.61, 67.12)
- alpha    ~ Gamma(3.267, 0.109)
"""
function sample_substitution_params(seed) 
    rng = StableRNG(seed)   # stable RNG object
    # kappa ~ LogNormal
    kappa = rand(rng, LogNormal(1.4215, 0.2798))
    # base frequencies ~ Dirichlet
    freq_dist = Dirichlet([66.59, 38.41, 38.61, 67.12])
    basefreqs = rand(rng, freq_dist)
    # gamma shape alpha ~ Gamma
    alpha = rand(rng, Gamma(3.267, 0.109))
    return (kappa=kappa, basefreqs=basefreqs, alpha=alpha)
end

"""
    sample_substitution_params_plant(seed)

Same as `sample_substitution_params` for the plant (angiosperm) tree, with
distributions fitted to the angiosperm loci (speciestree_angiosperm.jl, step 7):
- kappa    ~ LogNormal(1.3330, 0.2079), from the loci for which ModelFinder chose HKY
- basefreqs ~ Dirichlet(101.99, 69.27, 90.29, 100.32), from all loci with F
- alpha    ~ LogNormal(-0.3262, 0.7122), from all loci with G
"""
function sample_substitution_params_plant(seed)
    rng = StableRNG(seed)
    kappa = rand(rng, LogNormal(1.3330, 0.2079))
    basefreqs = rand(rng, Dirichlet([101.99, 69.27, 90.29, 100.32]))
    alpha = rand(rng, LogNormal(-0.3262, 0.7122))
    return (kappa=kappa, basefreqs=basefreqs, alpha=alpha)
end

"""
Parse the parameter name into the dictionary used by generate_master_seed
Input: params_string (String): 
    DUP<dup_rate>-LOS<loss_rate>-RV<ratevar>-N_ind<n_inds>-SF<SF>[-genelen<genelen>]
    e.g. "DUP0.01-LOS0.02-RVG-N_ind10-SF1.0-genelen1000"
Output: Dict with keys "dup_rate" (Float32), "LOSS_RATE" (Float32), 
    "ratevar" (String), "n_inds" (Int),
    "scaling_factor" (Float32) and "genelen" (Float32, 
    or nothing for old runs without genelen in the name)
Notes: the string is split at "-LOS", "-RV", "-N_ind", "-SF", "-genelen" 
    instead of using a regex, because the rates can be in scientific 
    notation (1.0e-5) which has a "-" inside
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
        scaling_factor_str = split(part5, "-genelen")[1]
        part6 = split(part5, "-genelen")[2]
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
Mapping from individuals to species, from the tip labels of the gene trees 
    (tips are Species_Accession, e.g. A_1)
Input:
    -gene_trees (Vector{HybridNetwork})
    -output_dir (String): where to write the mapping file
    -mode (String): "snaq" (default) writes id_to_species.csv
        with columns individual, species, for snaq!
        "astral" writes astral_mapping.txt with one "individual species" per line, 
        for ASTRAL -a
Output: Tuple (Dict individual => species, path of the mapping file)
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

Replace `:nan` in ASTRAL newick output with `:0.000001` so PhyloNetworks
can parse the string without errors.
"""
function clean_newick_nan(newick_str::String)
    # Replace :nan with :0.000001 (a small default value)
    return replace(newick_str, ":nan" => ":0.000001")
end

"""
    simplify_tip_labels(tree::HybridNetwork) -> HybridNetwork

Strip accession IDs from tip labels in-place, keeping only the prefix
before the first underscore (e.g. "A_0" -> "A"). Modifies the network.
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
# Retrieved from SNaQ.jl multiplealleles docs.
"""
    hasrep

Return true if a row (4-taxon set) has a "repeated" species, that is, a species
whose name ends with "__2". Otherwise, return false.

Warning: this function assumes that taxon names are in columns
"t1", "t2", "t3", "t4". For data frames with different column names,
e.g. "taxon1", "taxon2" etc., simply edit the code below by replacing
`:t1` by `:taxon1` 
"""
function hasrep(row)
    occursin(r"__2$", row[:t1]) || occursin(r"__2$", row[:t2]) ||
    occursin(r"__2$", row[:t3]) || occursin(r"__2$", row[:t4])
end 

"""
set_up_paramname_root(dup_rate, loss_rate, ratevar,
            n_inds, scaling_factor_branch_length,
            gene_len) -> String

Build the canonical parameter name used in directory and file naming.
Format: `DUP<d>-LOS<l>-RV<r>-N_ind<n>-SF<s>-genelen<g>`
"""
function set_up_paramname_root(dup_rate, loss_rate,
        ratevar, n_inds, scaling_factor_branch_length, gene_len)
    return "DUP$dup_rate-LOS$loss_rate-RV$ratevar" *
        "-N_ind$n_inds-SF$scaling_factor_branch_length-genelen$gene_len"
end


"""
create_mapping_file_based_on_gene_tree(mapping_file::String)
        -> Dict{String, String}
Read an astral mapping file (one "individual species" pair per line,
whitespace separated) into a dictionary individual => species.
"""
function create_mapping_file_based_on_gene_tree(mapping_file::String)
    individual_to_species = Dict{String, String}()
    for line in readlines(mapping_file)
        parts = split(strip(line))
        if length(parts) == 2
            individual_to_species[parts[1]] = parts[2]
        end
    end
    return individual_to_species
end

"""
    modify_locus_tree_labels!(locus_tree::HybridNetwork) -> HybridNetwork
Rename the tips of a locus tree in place, keeping only the part before the
first underscore (e.g. "A_0" -> "A"), so they match the species tree.
"""
function modify_locus_tree_labels!(locus_tree::HybridNetwork)
    for tip in locus_tree.leaf
        tip.name = String(split(tip.name, "_")[1])
    end
    return locus_tree
end

"""
map_gene_tree_based_on_mapping_file(
        gene_tree::HybridNetwork, 
        individual_to_species::Dict{String, String}
    ) -> HybridNetwork
Map gene tree tip labels from individuals (e.g. "A_1", "B_2") to species
names (e.g. "A", "B") using the provided dictionary. Returns a modified copy.
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
Group the tips of a tree by species: the species is the part of the tip name before the "_"
Input: tree (HybridNetwork)
Output: Dict species => vector of tip names
Example: dict = group_tips_by_species(net)
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
subtree_for_combination(tree::HybridNetwork, 
        tips_to_keep::Vector{String}) -> HybridNetwork 
Create a subtree by keeping only specified tips from the original tree.
This function takes a phylogenetic tree and a list of tip labels to retain.
It returns a new tree that includes only those tips, pruning all others.    
"""
function subtree_for_combination(tree, tips_to_keep)
    copy_tree = deepcopy(tree)
    tips_to_drop = setdiff(tipLabels(copy_tree), tips_to_keep)
    for t in tips_to_drop
        deleteleaf!(copy_tree, t) 
    end
    return copy_tree
end

"""
Length of the first sequence in a fasta alignment
    (all sequences have the same length in an alignment)
Input: fasta_file (String)
Output: Int, number of sites including gaps. 
    Sequences on several lines are fine, we stop after the first sequence
"""
function get_alignment_length(fasta_file::String)
    open(fasta_file, "r") do f
        current_seq = ""
        started = false
        
        for line in eachline(f)
            if startswith(line, ">")
                if started && !isempty(current_seq) # next header: done
                    return length(current_seq)
                end
                started = true
            elseif started
                current_seq *= strip(line)
            end
        end
        return length(current_seq)
    end
end 

"""
Why did simphy fail? Look at its l_trees.trees for missing file,
    trees with a single species,
or fewer than 4 species left after removing the "Lost-" tips
Input: simphy_output_dir (String), rep and iter (Int, for the message)
Output: String with the diagnosis, to print in the log
"""
function diagnose_simphy_locus_trees(simphy_output_dir::String,
        rep::Int, iter::Int)
    # simphy may write in a subfolder named "1"
    l_trees_path = joinpath(simphy_output_dir, "l_trees.trees")
    if !isfile(l_trees_path)
        l_trees_path_alt = joinpath(simphy_output_dir, "1", "l_trees.trees")
        if !isfile(l_trees_path_alt)
            return "SimPhy failed before generating locus tree file" *
                " (l_trees.trees not found)"
        end
        l_trees_path = l_trees_path_alt
    end
    
    locus_trees = readlines(l_trees_path)
    num_locus_trees = length(locus_trees)
    problematic_trees = [] # (locus id, reason)
    for (locus_id, newick_str) in enumerate(locus_trees)
        if isempty(strip(newick_str))
            push!(problematic_trees, (locus_id, "Empty tree string"))
            continue
        end
        
        try
            net = readnewick(newick_str)
            all_tips = [tip.name for tip in net.leaf]
            # remove the tips lost by simphy
            lost_tips = filter(name -> startswith(name, "Lost-"), all_tips)
            for lost_tip in lost_tips
                try
                    deleteleaf!(net, lost_tip)
                catch e
                    push!(problematic_trees,
                        (locus_id, "Failed to remove Lost tip: $lost_tip"))
                    continue
                end
            end
            
            remaining_tips = [tip.name for tip in net.leaf]
            num_remaining = length(remaining_tips)
            species_list = [] # species = part before the first "_"
            for tip_name in remaining_tips
                parts = split(tip_name, "_")
                species = parts[1]
                push!(species_list, species)
            end
            unique_species = unique(species_list)
            num_unique_species = length(unique_species)

            if num_remaining == 0
                push!(problematic_trees,
                    (locus_id, "All taxa were Lost (no remaining taxa)"))
            elseif num_remaining == 1
                species = unique_species[1]
                push!(problematic_trees,
                    (locus_id,
                    "Only 1 taxon ($species) remains after removing Lost taxa"))
            elseif num_unique_species == 1 # several tips, one species
                species = unique_species[1]
                push!(problematic_trees,
                    (locus_id,
                    "All $num_remaining remaining taxa are from species " *
                    "'$species' ($(join(remaining_tips, ", ")))"))
            elseif num_unique_species < 4
                push!(problematic_trees,
                    (locus_id,
                    "Only $num_unique_species species remain: " *
                    "$(join(unique_species, ", ")) " *
                    "(need >=4 for phylogenetic inference)"))
            end

        catch e
            push!(problematic_trees,
                (locus_id, "Failed to parse or process tree: $e"))
        end
    end

    if isempty(problematic_trees)
        return "SimPhy failed after generating $num_locus_trees locus trees," *
            " but cause is unknown. All locus trees appear valid." *
            " Inspect the locus tree file: $l_trees_path"
    else
        n_bad = length(problematic_trees)
        msg = "SimPhy generated $num_locus_trees locus trees" *
            " but $n_bad are problematic:\n"
        for (locus_id, reason) in problematic_trees
            msg *= "  Locus tree #$locus_id: $reason\n"
        end
        msg *= "Locus tree file: $l_trees_path"
        return msg
    end
end
