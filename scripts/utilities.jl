#= Utility functions shared across all pipeline scripts.

Part 0 -- Logging helpers for distributed workers
    worker_println  -> println that also writes to worker log file
    open_worker_log -> open per-worker log at start of parallel job
    close_worker_log -> close per-worker log at end of parallel job

Part 1 -- Gene tree cleanup after SimPhy
    replace_taxa_name      -> rename one tip given species char + newick
    count_missing          -> count missing/duplicate taxa in a newick
    modify_newicks         -> filter and simplify a gene tree file
    modify_newick_for_n_genes -> run modify_newicks across n gene files

Part 2 -- Simulation control helpers
    seed_generator   -> build N x M seed arrays from a master seed
    calculate_batch  -> estimate how many SimPhy iters remain

Part 3 -- General utilities
    pad_number             -> zero-padded string matching SimPhy output
    check_existing_dir     -> interactive prompt to clear output dirs
    replace_tips_with_letters -> convert species names to letters
    replace_pop_inInd_file -> assign unique pop IDs in .ind files
    generate_software_seeds -> per-software seeds from master seed
    generate_master_seed   -> deterministic seed from param dict
    get_dict_for_seed_setting -> extract seed params from name string
    totallength            -> sum all branch lengths in a tree/network
    hasrep                 -> detect repeated taxa in a newick string
    set_up_paramname_root  -> build the canonical param name string

Part 4 -- Summary statistics
    parse_parameter_setting      -> parse param name string to dict
    summarize_gamma_by_threshold -> gamma stats grouped by parameter
    summarize_WR_by_threshold    -> WR stats grouped by parameter
=#
using PhyloNetworks
using StatsBase
using StableRNGs
using Primes
using Statistics
using Printf
using InlineStrings
using IterTools
using Glob
using Distributions

#-----------------------------------------------#       
#               Part 0  
#    Logging helper functions
#-----------------------------------------------#
"""
    worker_println(args...)

Print to stdout and flush. If logging is enabled, also write to the
worker-specific log file.

Globals needed: `log_output` (Bool), `worker_logfile` (IOStream or nothing).
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
            # Silent fail if log write fails
        end
    end
end

"""
    open_worker_log()

Open `screen_<paramname>_worker<id>.log` in `outfolder`.
Call at the start of each distributed worker.

Globals needed: `log_output`, `outfolder`, `paramname_root`.
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
    close_worker_log()

Close the worker log file opened by `open_worker_log()`.
Call at the end of each distributed worker. Silent on failure.
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

Count and optionally rename tips for one species in a SimPhy newick.
Tips are labelled `char_locusID_accessionID`; this strips the locusID.

Returns `(count, modified_newick)`:
- count == 0: species absent, newick unchanged.
- count == 1: one unique copy, locus ID stripped.
- count >= 2: duplicated gene copy, newick unchanged.
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
    count_missing(newick::String)

Check a SimPhy newick (taxa A–H, format `char_x_y`) for missing or
duplicate gene copies. Returns `(n_missing, simplified_newick)` if all
copies are unique, or `nothing` if any taxon is duplicated (paralogy not
hidden).
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
    modify_newicks(input::String, output::String) -> Tuple{Int, Int}

Filter and simplify gene trees from a SimPhy output file.
Trees with duplicated gene copies or more than `max_taxa_missing`
absent taxa are dropped. Valid trees have tip labels simplified
from `A_locusID_accessionID` to `A_accessionID`.

Returns `(n_repeated, n_missing)` — counts of excluded trees.
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
    
    # Read atomically to avoid AFS race conditions where open/read/close
    # can interleave with concurrent writes. Each SimPhy file = one newick.
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
        @error "Failed to process tree in $input" exception=(e,
            catch_backtrace())
        return (0, 0)
    end

    # helper to remove a file if it exists, with a warning on failure
    function try_rm(path)
        isempty(path) && return
        try
            isfile(path) && rm(path; force=true)
        catch e
            @warn "Failed to remove $path" exception=(e, catch_backtrace())
        end
    end

    if isnothing(temp)
        # Tree has repeated taxa - invalid
        println("The input $input has repeated taxa")

        if !debug_mode
            try_rm(intermediate_file_mapsl)
            try_rm(intermediate_file_maplg)
            try_rm(input)
        else
            println("  DEBUG: Keeping intermediate files and tree for analysis")
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
            @warn "Failed to write $output" exception=(e, catch_backtrace())
        end
    end
    
    # Return counts of invalid trees
    return tree_repeated_taxa, tree_insufficient_taxa  
end


"""
    modify_newick_for_n_genes(input, output, n, iteration) -> Tuple{Int, Int}

Run `modify_newicks` across `n` gene tree files from one SimPhy replicate.
Input files: `g_trees{GENEID}.trees` (zero-padded IDs).
Output files: `g_trees_noLocusID_Gene{GENEID}_Int{iteration}.trees`.

Returns `(n_repeated, n_missing)` totalled across all genes.
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
        
        # search maplg only — when dup=0, only maplg exists (no mapsl)
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

        # Move output files one level up, removing the nested Int/1 folder
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
    seed_generator(master_seed, n, m, output_dir, output_file_name)
        -> Matrix{Int32}

Generate an n×m seed matrix via StableRNG and write it as a tab-separated
file (header: col1…colm, one rep per row). Returns the matrix.
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

Estimate how many gene trees to simulate next given `num` currently valid
trees and a target of `n_genes`. Accounts for expected missing-tree rate.
Capped at 3× n_genes to avoid excessive resource usage.

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
    pad_number(num::Int, range::Int) -> String

Zero-pad `num` to as many digits as `range` requires.
Used for consistent SimPhy-style filenames (e.g. gene 1 of 200 → "001").
"""
function pad_number(num::Int, range::Int) 
  digits = ceil(Int, log10(range + 1))
  number_string = lpad(string(num), digits, '0') 
  return number_string
end

"""
    get_hybrid_info(n::HybridNetwork, outgroup::String="A")

For an H=1 network, return the leaf sets at the hybrid node rooted on
`outgroup`. Result: `(hybrid_taxon, major_donor, minor_donor)` as
comma-joined strings (e.g. `"B,C"`). Returns "NA" fields on any error.

`n` is mutated in-place (rerooted). Pass `deepcopy(n)` to preserve original.
Falls back to major-edge rooting if `rootatnode!` raises `RootMismatch`.

# Examples
```julia-repl
julia> net = readnewick("(E,((((H,(G)#H1:::0.93)...);")
julia> info = get_hybrid_info(net, "A")
julia> info.hybrid_taxon   # e.g. "G"
julia> info.major_donor    # e.g. "B,C,D,E,F,H"
julia> info.minor_donor    # e.g. "A"
```
"""
function get_hybrid_info(n::HybridNetwork, outgroup::String="A")

    # Downward-only traversal: collect leaves below start_node,
    # never going through stop_node.
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

    # ----------------------------------------------------------------
    # KEY FIX: start from `start_node` and try leaves_below(., hnode).
    # If empty (because the donor parent's only child IS hnode), climb
    # up via tree edges until a node with reachable leaves is found.
    # This correctly scopes to the donor clade without spilling into
    # the entire network (which the undirected traversal did).
    #
    # Example: 
    #   admixb --only child--> #H1  =>  climb to rooty
    #   leaves_below(rooty, #H1)  = {H}           [major donor] v
    #   admixa --only child--> #H1  =>  climb to G_z
    #   leaves_below(G_z,    #H1)  = {B, C}        [minor donor] v
    # An example when leabes_below would fail: 
    # (((F:0.144,((D:0.077,E:0.075)Rrrl:0.045,((B:0.04,C:0.04)B_a:0.054,
    # (#H1:0.0::0.054)admixa:0.0)G_z:0.004)Rrrz:0.012)Rrr:0.021,
    # (H:0.159,((G:0.081)#H1:0.0::0.946)admixb:0.072)rooty:0.134)
    # root:0.126,A:0.126)rootb;
    # ----------------------------------------------------------------
    function donor_leaves(start_node, hnode)
        current = start_node
        while true
            lv = leaves_below(current, hnode)
            isempty(lv) || return lv

            # Climb via the tree-parent edge only (skip hybrid edges
            # to avoid following a reticulation upward).
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
        # --- Root on outgroup ---
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

        # --- Validate single hybrid ---
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
    replace_tips_with_letters(tree::Union{HybridNetwork, String})
        -> Union{HybridNetwork, String}

Replace tip labels with capital letters A–Z (by order of appearance).
Handles both Newick strings and HybridNetwork objects.
Returns original with a warning if there are more than 26 tips.

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
    replace_pop_inInd_file(input_file::String, output_file::String)

Replace the "POP" placeholder in column 3 of an eigenstrat `.ind` file
with the species letter extracted from column 1 (e.g. "A_0" → "A").
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
    - `"genelen"`: `Float32` or `nothing` — gene length, if present

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
    `"individual species"` (`astral_mapping.txt`), for ASTRAL’s `-a` flag.

# Returns
- `Tuple{Dict{String,String}, String}`: A tuple with:
  1. A dictionary mapping each accession (e.g., `"A_1"`) to its species
     (e.g., `"A"`).
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
before the first underscore (e.g. "A_0" → "A"). Modifies the network.
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
Iterates over a set of examples.

# Description
This loop is intended to process or analyze each example in a collection. 

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
        deleteleaf!(copy_tree, t) 
    end
    return copy_tree
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
    diagnose_simphy_locus_trees(simphy_output_dir, rep, iter) -> String

Scan SimPhy's `l_trees.trees` and return a diagnostic string explaining
why gene tree generation failed. Checks for missing file, single-species
trees, and insufficient diversity after removing "Lost-" tips.
"""
function diagnose_simphy_locus_trees(simphy_output_dir::String,
        rep::Int, iter::Int)
    # Check if l_trees.trees file exists
    # SimPhy may create output in a subdirectory named "1"
    l_trees_path = joinpath(simphy_output_dir, "l_trees.trees")
    if !isfile(l_trees_path)
        # Try the "1" subdirectory
        l_trees_path_alt = joinpath(simphy_output_dir, "1", "l_trees.trees")
        if !isfile(l_trees_path_alt)
            return "SimPhy failed before generating locus tree file" *
                " (l_trees.trees not found)"
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
                    push!(problematic_trees,
                        (locus_id, "Failed to remove Lost tip: $lost_tip"))
                    continue
                end
            end
            
            # Get remaining tips after removing "Lost-" taxa
            remaining_tips = [tip.name for tip in net.leaf]
            num_remaining = length(remaining_tips)
            
            # Extract species prefix (before first underscore)
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
                push!(problematic_trees,
                    (locus_id, "All taxa were Lost (no remaining taxa)"))
            elseif num_remaining == 1
                # Only one taxon remains
                species = unique_species[1]
                push!(problematic_trees,
                    (locus_id,
                    "Only 1 taxon ($species) remains after removing Lost taxa"))
            elseif num_unique_species == 1
                # Multiple taxa but all from same species
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
                    "(need ≥4 for phylogenetic inference)"))
            end

        catch e
            push!(problematic_trees,
                (locus_id, "Failed to parse or process tree: $e"))
        end
    end

    # Generate diagnostic message
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


# ============================================================
# Part 4 -- Summary statistics
# ============================================================

"""
    parse_parameter_setting(parameter_setting::String)

Parse parameter setting string to extract key parameters.
Returns a dictionary with keys: ratevar, duploss_rate, SF, n_inds

# Example
- Input: "DUP0.0003-LOS0.0003-RVN-N_ind2-SF1.0-genelen1000"
- Output: `Dict("ratevar" => "RVN", "duploss_rate" => "0.0003",
           "SF" => "1.0", "n_inds" => "2")`
"""
function parse_parameter_setting(parameter_setting::String)
    params = Dict{String, String}()

    # Extract ratevar
    m = match(r"-(RVG|RVL|RVN|RVGL)-", parameter_setting)
    params["ratevar"] = m !== nothing ? m.captures[1] : "unknown"

    # Extract duplication/loss rate
    m = match(r"DUP([\d\.e\-]+)-LOS([\d\.e\-]+)", parameter_setting)
    params["duploss_rate"] = m !== nothing ? m.captures[1] : "unknown"

    # Extract SF (scaling factor)
    m = match(r"SF([\d\.]+)", parameter_setting)
    params["SF"] = m !== nothing ? m.captures[1] : "unknown"

    # Extract n_inds (number of individuals)
    m = match(r"N_ind(\d+)", parameter_setting)
    params["n_inds"] = m !== nothing ? m.captures[1] : "unknown"

    return params
end


"""
    summarize_gamma_by_threshold(input_dir::String,
                                 column_name_1::Symbol,
                                 column_name_2::Symbol,
                                 output_dir::String,
                                 output_filename::String,
                                 thresholds::Vector{Float64}=[0.05, 0.1, 0.25])

Summary statistics for gamma values below given thresholds,
grouped by parameter.

# Arguments
- `input_dir`: Directory containing CSV files
- `column_name_1`: First gamma column
- `column_name_2`: Second gamma column
- `output_dir`: Directory to save output file
- `output_filename`: Name of output file
- `thresholds`: Vector of threshold values (default: [0.05, 0.1, 0.25])
"""
function summarize_gamma_by_threshold(input_dir::String,
        column_name_1::Symbol, column_name_2::Symbol,
        output_dir::String, output_filename::String,
        thresholds::Vector{Float64}=[0.05, 0.1, 0.25])
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))

    if isempty(csv_files)
        println("No CSV files found in $input_dir")
        return
    end

    mkpath(output_dir)

    # Store data for each category
    category_stats = Dict{String, Dict{String, Any}}()
    categories = ["ratevar", "duploss_rate", "SF", "n_inds"]

    for category in categories
        category_stats[category] = Dict{String, Vector{Float64}}()
    end

    # Overall statistics
    all_gamma_values = Float64[]

    # Process each CSV file
    for csv_file in csv_files
        filepath = joinpath(input_dir, csv_file)
        df = CSV.read(filepath, DataFrame)

        if !hasproperty(df, column_name_1) || !hasproperty(df, column_name_2)
            println("Warning: Columns not found in $csv_file")
            continue
        end

        parameter_setting = replace(csv_file, r"^(SNaQ-|findgraph-)" => "")
        parameter_setting = replace(
            parameter_setting, r"(-summary)?\.csv$" => "")

        # Parse parameters
        params = parse_parameter_setting(parameter_setting)

        # Calculate minor gamma
        minor_gamma = min.(df[!, column_name_1], df[!, column_name_2])
        minor_gamma = filter(!isnan, collect(skipmissing(minor_gamma)))

        if isempty(minor_gamma)
            continue
        end

        # Add to overall statistics
        append!(all_gamma_values, minor_gamma)

        # Add to category-specific statistics
        for category in categories
            key = params[category]
            if haskey(category_stats[category], key)
                append!(category_stats[category][key], minor_gamma)
            else
                category_stats[category][key] = copy(minor_gamma)
            end
        end
    end

    # Write summary to file
    output_path = joinpath(output_dir, output_filename)
    open(output_path, "w") do io
        println(io, "Gamma summary statistics")
        println(io, "Thresholds: ", join(thresholds, ", "))
        println(io, "="^70)
        println(io)

        # Overall statistics
        n_total = length(all_gamma_values)

        println(io, "Overall statistics")
        println(io, "-"^70)
        println(io, "Total gamma values: $n_total")
        println(io, "Mean: $(round(mean(all_gamma_values), digits=3))")
        println(io, "Median: $(round(median(all_gamma_values), digits=3))")
        println(io, "Std dev: $(round(std(all_gamma_values), digits=3))")
        println(io)

        for threshold in thresholds
            n_below = count(x -> x < threshold, all_gamma_values)
            pct_below = n_total > 0 ?
                round(n_below / n_total * 100, digits=1) : 0.0
            println(io, "  Below $threshold: $n_below ($pct_below%)")
        end
        println(io)
        println(io)

        # Category-specific statistics
        for category in categories
            println(io, "Statistics by $category")
            println(io, "-"^70)

            category_data = category_stats[category]
            sorted_keys = sort(collect(keys(category_data)))

            for key in sorted_keys
                data = category_data[key]
                n_total_cat = length(data)
                mean_gamma = round(mean(data), digits=3)
                median_gamma = round(median(data), digits=3)
                std_gamma = round(std(data), digits=3)

                println(io, "  $key:")
                println(io,
                    "    n=$n_total_cat, mean=$mean_gamma," *
                    " median=$median_gamma, std=$std_gamma")

                for threshold in thresholds
                    n_below_cat = count(x -> x < threshold, data)
                    pct_below_cat = n_total_cat > 0 ?
                        round(n_below_cat / n_total_cat * 100, digits=2) : 0.0
                    println(io,
                        "    <$threshold: $n_below_cat ($pct_below_cat%)")
                end
                println(io)
            end

            println(io)
        end

        println(io, "="^70)
        println(io, "Source: $input_dir")
    end

    println("Gamma summary saved to: $output_path")

    # Create CSV table with all statistics
    csv_rows = []

    # Add overall statistics row
    n_total = length(all_gamma_values)
    mean_gamma = round(mean(all_gamma_values), digits=4)
    median_gamma = round(median(all_gamma_values), digits=4)
    sd_gamma = round(std(all_gamma_values), digits=4)

    overall_row = Dict(
        "category" => "Overall",
        "value" => "all",
        "n" => n_total,
        "mean_gamma" => mean_gamma,
        "median_gamma" => median_gamma,
        "sd_gamma" => sd_gamma
    )

    for threshold in thresholds
        n_below = count(x -> x < threshold, all_gamma_values)
        pct_below = n_total > 0 ?
            round(n_below / n_total * 100, digits=1) : 0.0
        col_name = "gamma_below_$(threshold)"
        overall_row[col_name] = pct_below
    end

    push!(csv_rows, overall_row)

    # Add category-specific statistics
    for category in categories
        category_data = category_stats[category]
        sorted_keys = sort(collect(keys(category_data)))

        for key in sorted_keys
            data = category_data[key]
            n_total_cat = length(data)
            mean_gamma_cat = round(mean(data), digits=3)
            median_gamma_cat = round(median(data), digits=3)
            sd_gamma_cat = round(std(data), digits=3)

            row = Dict(
                "category" => category,
                "value" => key,
                "n" => n_total_cat,
                "mean_gamma" => mean_gamma_cat,
                "median_gamma" => median_gamma_cat,
                "sd_gamma" => sd_gamma_cat
            )

            for threshold in thresholds
                n_below_cat = count(x -> x < threshold, data)
                pct_below_cat = n_total_cat > 0 ?
                    round(n_below_cat / n_total_cat * 100, digits=1) : 0.0
                col_name = "gamma_below_$(threshold)"
                row[col_name] = pct_below_cat
            end

            push!(csv_rows, row)
        end
    end

    # Convert to DataFrame and save as CSV
    summary_df = DataFrame(csv_rows)

    # Reorder columns for better readability
    col_order = [
        "category", "value", "n",
        "mean_gamma", "median_gamma", "sd_gamma"]
    for threshold in thresholds
        push!(col_order, "gamma_below_$(threshold)")
    end
    select!(summary_df, col_order)

    # Generate CSV filename
    csv_filename = replace(output_filename, r"\.(txt|log)$" => ".csv")
    csv_output_path = joinpath(output_dir, csv_filename)

    CSV.write(csv_output_path, summary_df)
    println("Gamma summary CSV saved to: $csv_output_path")
end


"""
    summarize_WR_by_threshold(input_dir::String,
                              column_name::Symbol,
                              output_dir::String,
                              output_filename::String,
                              metric_name::String,
                              percentiles::Vector{Float64}=[0.95, 0.99])

Summary statistics for worst residual (WR) values grouped by parameters.
Works for single-column metrics like H0_best_tree_WR and H1_best_graph_WR.

# Arguments
- `input_dir`: Directory containing CSV files
- `column_name`: Column to analyze (e.g., :H0_best_tree_WR)
- `output_dir`: Directory to save output file
- `output_filename`: Name of output file
- `metric_name`: Display name for the metric (e.g., "H0 Best Tree WR")
- `percentiles`: Percentiles to report (default: [0.95, 0.99])
"""
function summarize_WR_by_threshold(input_dir::String,
                                   column_name::Symbol,
                                   output_dir::String,
                                   output_filename::String,
                                   metric_name::String,
                                   percentiles::Vector{Float64}=[0.95, 0.99])
    csv_files = filter(x -> endswith(x, ".csv"), readdir(input_dir))

    if isempty(csv_files)
        println("No CSV files found in $input_dir")
        return
    end

    mkpath(output_dir)

    # Store data for each category
    category_stats = Dict{String, Dict{String, Any}}()
    categories = ["ratevar", "duploss_rate", "SF", "n_inds"]

    for category in categories
        category_stats[category] = Dict{String, Vector{Float64}}()
    end

    # Per-file percentile tracking (for min/max across parameter settings)
    # Structure: category -> key -> pct -> [per_file_pct_value, ...]
    category_per_file_pcts =
        Dict{String, Dict{String, Dict{Float64, Vector{Float64}}}}()
    for category in categories
        category_per_file_pcts[category] =
            Dict{String, Dict{Float64, Vector{Float64}}}()
    end
    all_per_file_pcts =
        Dict{Float64, Vector{Float64}}(pct => Float64[] for pct in percentiles)

    # ratevar x SF interaction
    ratevar_x_SF_stats = Dict{String, Vector{Float64}}()
    ratevar_x_SF_per_file_pcts =
        Dict{String, Dict{Float64, Vector{Float64}}}()

    # ratevar_group x ILS: (RVG+RVN vs RVL) × ILS level (from SF)
    ils_map = Dict("1.0" => "high", "0.5" => "low")
    ratevar_group_x_ILS_stats = Dict{String, Vector{Float64}}()
    ratevar_group_x_ILS_per_file_pcts =
        Dict{String, Dict{Float64, Vector{Float64}}}()

    # Overall statistics
    all_wr_values = Float64[]

    # Process each CSV file
    for csv_file in csv_files
        filepath = joinpath(input_dir, csv_file)
        df = CSV.read(filepath, DataFrame)

        if !hasproperty(df, column_name)
            println("Warning: Column $(column_name) not found in $csv_file")
            continue
        end

        parameter_setting = replace(csv_file, r"^(SNaQ-|findgraph-)" => "")
        parameter_setting = replace(
            parameter_setting, r"(-summary)?\.csv$" => "")

        # Parse parameters
        params = parse_parameter_setting(parameter_setting)

        # Extract WR values and filter NaN/missing
        wr_values = filter(!isnan, collect(skipmissing(df[!, column_name])))

        if isempty(wr_values)
            continue
        end

        # Add to overall statistics
        append!(all_wr_values, wr_values)

        # Compute per-file percentiles and accumulate for min/max tracking
        file_pcts =
            Dict(pct => quantile(wr_values, pct) for pct in percentiles)
        for pct in percentiles
            push!(all_per_file_pcts[pct], file_pcts[pct])
        end

        # Add to category-specific statistics
        for category in categories
            key = params[category]
            if haskey(category_stats[category], key)
                append!(category_stats[category][key], wr_values)
            else
                category_stats[category][key] = copy(wr_values)
            end
            # Per-file percentile tracking
            if !haskey(category_per_file_pcts[category], key)
                category_per_file_pcts[category][key] =
                    Dict(pct => Float64[] for pct in percentiles)
            end
            for pct in percentiles
                push!(category_per_file_pcts[category][key][pct],
                    file_pcts[pct])
            end
        end

        # ratevar x SF interaction
        rv_sf_key = "$(params["ratevar"])_SF$(params["SF"])"
        if haskey(ratevar_x_SF_stats, rv_sf_key)
            append!(ratevar_x_SF_stats[rv_sf_key], wr_values)
        else
            ratevar_x_SF_stats[rv_sf_key] = copy(wr_values)
        end
        if !haskey(ratevar_x_SF_per_file_pcts, rv_sf_key)
            ratevar_x_SF_per_file_pcts[rv_sf_key] =
                Dict(pct => Float64[] for pct in percentiles)
        end
        for pct in percentiles
            push!(ratevar_x_SF_per_file_pcts[rv_sf_key][pct], file_pcts[pct])
        end

        # ratevar_group x ILS: RVG+RVN vs RVL, crossed with ILS (SF)
        rv = params["ratevar"]
        rv_group = (rv == "RVL") ? "RVL" : "RVG+RVN"
        ils_level = get(ils_map, params["SF"], params["SF"])
        rv_ils_key = "$(rv_group)_ILS$(ils_level)"
        if haskey(ratevar_group_x_ILS_stats, rv_ils_key)
            append!(ratevar_group_x_ILS_stats[rv_ils_key], wr_values)
        else
            ratevar_group_x_ILS_stats[rv_ils_key] = copy(wr_values)
        end
        if !haskey(ratevar_group_x_ILS_per_file_pcts, rv_ils_key)
            ratevar_group_x_ILS_per_file_pcts[rv_ils_key] =
                Dict(pct => Float64[] for pct in percentiles)
        end
        for pct in percentiles
            push!(ratevar_group_x_ILS_per_file_pcts[rv_ils_key][pct],
                file_pcts[pct])
        end
    end

    # Write summary to file
    output_path = joinpath(output_dir, output_filename)
    open(output_path, "w") do io
        println(io, "$metric_name summary statistics")
        pct_labels = join(
            map(p -> "$(round(p*100, digits=0))%", percentiles), ", ")
        println(io, "Percentiles: ", pct_labels)
        println(io, "="^70)
        println(io)

        # Overall statistics
        n_total = length(all_wr_values)

        println(io, "Overall statistics")
        println(io, "-"^70)
        println(io, "Total WR values: $n_total")
        if !isempty(all_wr_values)
            println(io, "Mean: $(round(mean(all_wr_values), digits=4))")
            println(io, "Median: $(round(median(all_wr_values), digits=4))")
            println(io, "Std dev: $(round(std(all_wr_values), digits=4))")
        else
            println(io, "No valid data found for this metric")
        end
        println(io)

        # Calculate percentile values for overall data
        for pct in percentiles
            if !isempty(all_wr_values)
                q_value = round(quantile(all_wr_values, pct), digits=4)
                pct_label = round(Int, pct * 100)
                println(io,
                    "  $pct_label-th percentile: $q_value" *
                    " ($(pct_label)% of values below this)")
            end
        end
        println(io)

        # Per-setting percentile range (min/max across parameter settings)
        if !isempty(all_per_file_pcts[first(percentiles)])
            n_settings = length(all_per_file_pcts[first(percentiles)])
            println(io,
                "  Per-setting percentile range (min/max across" *
                " $n_settings parameter settings):")
            for pct in percentiles
                vals = all_per_file_pcts[pct]
                pct_label = round(Int, pct * 100)
                println(io,
                    "    $pct_label-th percentile:" *
                    " min=$(round(minimum(vals), digits=4))," *
                    " max=$(round(maximum(vals), digits=4))")
            end
        end
        println(io)
        println(io)

        # Category-specific statistics
        for category in categories
            println(io, "Statistics by $category")
            println(io, "-"^70)

            category_data = category_stats[category]
            sorted_keys = sort(collect(keys(category_data)))

            for key in sorted_keys
                data = category_data[key]
                n_total_cat = length(data)
                if !isempty(data)
                    mean_wr = round(mean(data), digits=4)
                    median_wr = round(median(data), digits=4)
                    std_wr = round(std(data), digits=4)

                    println(io, "  $key:")
                    println(io,
                        "    n=$n_total_cat, mean=$mean_wr," *
                        " median=$median_wr, std=$std_wr")

                    for pct in percentiles
                        if !isempty(data)
                            q_value = round(quantile(data, pct), digits=4)
                            pct_label = round(Int, pct * 100)
                            println(io,
                                "    $pct_label-th percentile (pooled):" *
                                " $q_value")
                        end
                    end

                    if haskey(category_per_file_pcts[category], key)
                        per_file = category_per_file_pcts[category][key]
                        n_settings = length(per_file[first(percentiles)])
                        println(io,
                            "    Per-setting range ($n_settings settings):")
                        for pct in percentiles
                            vals = per_file[pct]
                            pct_label = round(Int, pct * 100)
                            println(io,
                                "      $pct_label-th:" *
                                " min=$(round(minimum(vals), digits=4))," *
                                " max=$(round(maximum(vals), digits=4))")
                        end
                    end
                end
                println(io)
            end

            println(io)
        end

        # ratevar x SF interaction statistics
        println(io, "Statistics by ratevar x SF")
        println(io, "-"^70)

        for key in sort(collect(keys(ratevar_x_SF_stats)))
            data = ratevar_x_SF_stats[key]
            n_total_int = length(data)
            if !isempty(data)
                mean_wr = round(mean(data), digits=4)
                median_wr = round(median(data), digits=4)
                std_wr = round(std(data), digits=4)

                println(io, "  $key:")
                println(io,
                    "    n=$n_total_int, mean=$mean_wr," *
                    " median=$median_wr, std=$std_wr")

                for pct in percentiles
                    q_value = round(quantile(data, pct), digits=4)
                    pct_label = round(Int, pct * 100)
                    println(io,
                        "    $pct_label-th percentile (pooled): $q_value")
                end

                if haskey(ratevar_x_SF_per_file_pcts, key)
                    per_file = ratevar_x_SF_per_file_pcts[key]
                    n_settings = length(per_file[first(percentiles)])
                    println(io,
                        "    Per-setting range ($n_settings settings):")
                    for pct in percentiles
                        vals = per_file[pct]
                        pct_label = round(Int, pct * 100)
                        println(io,
                            "      $pct_label-th:" *
                            " min=$(round(minimum(vals), digits=4))," *
                            " max=$(round(maximum(vals), digits=4))")
                    end
                end
            end
            println(io)
        end

        # ratevar_group x ILS: (RVG+RVN vs RVL) × ILS level
        println(io, "Statistics by ratevar_group x ILS")
        println(io, "-"^70)

        for key in sort(collect(keys(ratevar_group_x_ILS_stats)))
            data = ratevar_group_x_ILS_stats[key]
            n_total_int = length(data)
            if !isempty(data)
                mean_wr = round(mean(data), digits=4)
                median_wr = round(median(data), digits=4)
                std_wr = round(std(data), digits=4)

                println(io, "  $key:")
                println(io,
                    "    n=$n_total_int, mean=$mean_wr," *
                    " median=$median_wr, std=$std_wr")

                for pct in percentiles
                    q_value = round(quantile(data, pct), digits=4)
                    pct_label = round(Int, pct * 100)
                    println(io,
                        "    $pct_label-th percentile (pooled): $q_value")
                end

                if haskey(ratevar_group_x_ILS_per_file_pcts, key)
                    per_file = ratevar_group_x_ILS_per_file_pcts[key]
                    n_settings = length(per_file[first(percentiles)])
                    println(io,
                        "    Per-setting range ($n_settings settings):")
                    for pct in percentiles
                        vals = per_file[pct]
                        pct_label = round(Int, pct * 100)
                        println(io,
                            "      $pct_label-th:" *
                            " min=$(round(minimum(vals), digits=4))," *
                            " max=$(round(maximum(vals), digits=4))")
                    end
                end
            end
            println(io)
        end

        println(io, "")
        println(io, "="^70)
        println(io, "Source: $input_dir")
    end

    println("$metric_name summary saved to: $output_path")

    # Create CSV table with all statistics
    csv_rows = []

    # Add overall statistics row
    n_total = length(all_wr_values)
    overall_row = Dict(
        "category" => "Overall",
        "value" => "all",
        "n" => n_total
    )

    if !isempty(all_wr_values)
        mean_wr = round(mean(all_wr_values), digits=4)
        median_wr = round(median(all_wr_values), digits=4)
        sd_wr = round(std(all_wr_values), digits=4)
        overall_row["mean_WR"] = mean_wr
        overall_row["median_WR"] = median_wr
        overall_row["sd_WR"] = sd_wr

        # Calculate percentile values for overall data
        for pct in percentiles
            q_value = round(quantile(all_wr_values, pct), digits=4)
            pct_label = round(Int, pct * 100)
            col_name = "percentile_$(pct_label)"
            overall_row[col_name] = q_value
        end

        # Per-setting min/max percentiles across all settings
        for pct in percentiles
            vals = all_per_file_pcts[pct]
            pct_label = round(Int, pct * 100)
            if !isempty(vals)
                overall_row["min_percentile_$(pct_label)"] =
                    round(minimum(vals), digits=4)
                overall_row["max_percentile_$(pct_label)"] =
                    round(maximum(vals), digits=4)
            end
        end
    end

    push!(csv_rows, overall_row)

    # Add category-specific statistics
    for category in categories
        category_data = category_stats[category]
        sorted_keys = sort(collect(keys(category_data)))

        for key in sorted_keys
            data = category_data[key]
            n_total_cat = length(data)

            row = Dict(
                "category" => category,
                "value" => key,
                "n" => n_total_cat
            )

            if !isempty(data)
                mean_wr_cat = round(mean(data), digits=4)
                median_wr_cat = round(median(data), digits=4)
                sd_wr_cat = round(std(data), digits=4)

                row["mean_WR"] = mean_wr_cat
                row["median_WR"] = median_wr_cat
                row["sd_WR"] = sd_wr_cat

                # Calculate pooled percentile values for this category
                for pct in percentiles
                    q_value = round(quantile(data, pct), digits=4)
                    pct_label = round(Int, pct * 100)
                    col_name = "percentile_$(pct_label)"
                    row[col_name] = q_value
                end

                # Per-setting min/max percentiles for this category value
                if haskey(category_per_file_pcts[category], key)
                    per_file = category_per_file_pcts[category][key]
                    for pct in percentiles
                        vals = per_file[pct]
                        pct_label = round(Int, pct * 100)
                        if !isempty(vals)
                            row["min_percentile_$(pct_label)"] =
                                round(minimum(vals), digits=4)
                            row["max_percentile_$(pct_label)"] =
                                round(maximum(vals), digits=4)
                        end
                    end
                end
            end

            push!(csv_rows, row)
        end
    end

    # Add ratevar x SF interaction rows
    for key in sort(collect(keys(ratevar_x_SF_stats)))
        data = ratevar_x_SF_stats[key]
        n_total_int = length(data)

        row = Dict(
            "category" => "ratevar_x_SF",
            "value" => key,
            "n" => n_total_int
        )

        if !isempty(data)
            row["mean_WR"]   = round(mean(data), digits=4)
            row["median_WR"] = round(median(data), digits=4)
            row["sd_WR"]     = round(std(data), digits=4)

            for pct in percentiles
                q_value   = round(quantile(data, pct), digits=4)
                pct_label = round(Int, pct * 100)
                row["percentile_$(pct_label)"] = q_value
            end

            if haskey(ratevar_x_SF_per_file_pcts, key)
                per_file = ratevar_x_SF_per_file_pcts[key]
                for pct in percentiles
                    vals      = per_file[pct]
                    pct_label = round(Int, pct * 100)
                    if !isempty(vals)
                        row["min_percentile_$(pct_label)"] =
                            round(minimum(vals), digits=4)
                        row["max_percentile_$(pct_label)"] =
                            round(maximum(vals), digits=4)
                    end
                end
            end
        end

        push!(csv_rows, row)
    end

    # Add ratevar_group x ILS interaction rows
    for key in sort(collect(keys(ratevar_group_x_ILS_stats)))
        data = ratevar_group_x_ILS_stats[key]
        n_total_int = length(data)

        row = Dict(
            "category" => "ratevar_group_x_ILS",
            "value" => key,
            "n" => n_total_int
        )

        if !isempty(data)
            row["mean_WR"]   = round(mean(data), digits=4)
            row["median_WR"] = round(median(data), digits=4)
            row["sd_WR"]     = round(std(data), digits=4)

            for pct in percentiles
                q_value   = round(quantile(data, pct), digits=4)
                pct_label = round(Int, pct * 100)
                row["percentile_$(pct_label)"] = q_value
            end

            if haskey(ratevar_group_x_ILS_per_file_pcts, key)
                per_file = ratevar_group_x_ILS_per_file_pcts[key]
                for pct in percentiles
                    vals      = per_file[pct]
                    pct_label = round(Int, pct * 100)
                    if !isempty(vals)
                        row["min_percentile_$(pct_label)"] =
                            round(minimum(vals), digits=4)
                        row["max_percentile_$(pct_label)"] =
                            round(maximum(vals), digits=4)
                    end
                end
            end
        end

        push!(csv_rows, row)
    end

    # Convert to DataFrame and save as CSV
    summary_df = DataFrame(csv_rows)

    # Reorder columns for better readability
    col_order = ["category", "value", "n"]
    if hasproperty(summary_df, :mean_WR)
        push!(col_order, "mean_WR")
    end
    if hasproperty(summary_df, :median_WR)
        push!(col_order, "median_WR")
    end
    if hasproperty(summary_df, :sd_WR)
        push!(col_order, "sd_WR")
    end
    for pct in percentiles
        pct_label = round(Int, pct * 100)
        col_name = "percentile_$(pct_label)"
        if hasproperty(summary_df, Symbol(col_name))
            push!(col_order, col_name)
        end
    end
    for pct in percentiles
        pct_label = round(Int, pct * 100)
        for prefix in ("min_percentile_", "max_percentile_")
            col_name = "$(prefix)$(pct_label)"
            if hasproperty(summary_df, Symbol(col_name))
                push!(col_order, col_name)
            end
        end
    end
    if length(col_order) > 3
        select!(summary_df, col_order)
    end

    # Generate CSV filename
    csv_filename = replace(output_filename, r"\.(txt|log)$" => ".csv")
    csv_output_path = joinpath(output_dir, csv_filename)

    CSV.write(csv_output_path, summary_df)
    println("$metric_name CSV summary saved to: $csv_output_path")
end
