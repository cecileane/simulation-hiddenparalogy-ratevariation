#= Script to store utility functions 
The below script has three parts 

Part 1 --> A group of unctions to modify gene tree output from SimPhy (scripts/simulation.jl)
    F1: "replace_taxa_name" -> Count and modify one tree tip given Char and newick 
    F2: "count_missing" -> Count # missing taxa in the given newick and modify all tree tips (conditions seen below)
    F3: "modify_newicks" -> write modified newicks (condtions seen below)
    F4: "modify_newick_for_n_genes" -> run "modify_newicks" for n times and output the files with given prefix 

Part 2 --> A group of helper functions to re-run simphy simulation and modify newicks (seen details in scripts/speciestree.jl section "run SimPhy + modify newick strings" ) 
    F1: “seed_generator” -> generate N x M arrays with random seeds given a master seed 
    F2: "calculate_batch" -> Estimate batch size for simulating gene trees in the next iterations 
    F3: "pad_number" -> convert a given Int to a string to meet the output from Simphy

Part 3 --> Other utilities functions 
    F1: "check_existing_dir" -> check if a path exists and then decide to remove the path or not
    F2: "replace_tips_with_letters" -> Replace text tips with letters to both strings and Network/tree (used in scripts/speciestree.jl) 
    F3: "replace_pop_inInd_files" -> Replace "POP" from the output ind file to make each sample is mapped to a unique population ID
=# 
using PhyloNetworks # Used in "replace_tips_with_letters" to get trees and networks 
using StatsBase # Use the countmap function to find duplicated items
using StableRNGs # Used to generate stanble random seeds that won't change between julia versions
using Primes # use "nextprimes()" function to get the next prime bigger than a number 

#-----------------------------------------------#       
#               Part 1  
#    Modify newick strings generated from Simphy
#    Goal: Mimic hidden paralogy    
#-----------------------------------------------#
"""
Inputs:
    -char (String): The target character to search for in tip labels.
    -tree_str (String): The Newick tree string containing tip labels in the format char_x_y.
Outputs: Tuple (count, modified_tree):
    -count (Integer): Number of occurrences of char in the tree.
    -modified_tree (String): The Newick string, modified if applicable.
    1) If tip appears in the string for only once --> Return (count, modified tree).
    2) If tip appears in the string for 0 or more than once --> Return (count, unmodified tree). 
Examples: 
# Case 1: Char appears in the string for only once
Input: replace_taxa_name('A', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
Output: (1, "(A_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")

# Case 2: Char appears in the string for 0 time
Input: replace_taxa_name('F', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
Output: (0, "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")

# Case 3: Char appears in the string for >= 2 times: 
Input: replace_taxa_name('A', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_1_0:0.02,A_2_1:0.02):0.01)));")
Output: (2, "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_1_0:0.02,A_2_1:0.02):0.01)));") # Note there are A_0_1, A_1_0, A_2_1. It output 2 here because A_0_1 and A_2_1 are gene copies from the same accesions. 

# Case 4: Similar to above, char appears in the string for >= 2 times, but the there are different individuals from the same species: 
Input: replace_taxa_name('A', "(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_0_0:0.02,A_0_2:0.02):0.01)));")
Output： (1, "(A_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(A_0:0.02,A_2:0.02):0.01)));")
becase all three As are the same gene copies from different accessions. 
"""
function replace_taxa_name(taxon::Char, newick::String)
    tre = readTopology(newick)
    tips = tipLabels(tre)

    # replace "A_(\d+)_(\d+)" with "A_\2"
    # From simphy tutorial, tips are labelled as taxa_locusID_accessionID, here, remove the locusID while keeping the accessionID. 
    regex = Regex("($(taxon))_\\d+_(\\d+)")
    matching_tips = filter(tip -> occursin(regex, tip), tips)

    # Needs to handle if there are multiple accesions in one species 
    # eg 1: A_1 and A_0 are different accessions for one species, so keep both 
    # eg 2: A_0 and A_0 are the same accession for one species that the genes got duplicated
    tip_list = [replace(tip, regex => s"\1_\2") for tip in matching_tips]
    counts = values(countmap(tip_list)) 
    num = isempty(counts) ? 0 : maximum(counts) # Only take the maximum. Thus, if one individual has more than one copies, the tree will be discarded later. If the char is not in the string, num = 0 

    if num == 0 || num >= 2 # If # of taxons >= 2 and == 0, no change to the tips: 
        return (num, newick)
    else # If # of taxons == 1, remove locus ID from the tip 
        modified_newick = replace(newick, regex => s"\1_\2") 
        return(num, modified_newick)
    end
end

"""
Inputs:
    -tree_str (String): A Newick tree string containing taxa in the format char_x_y. Taxa are hardcoded as "ABCDEFGH".
Outputs:
    -Tuple (missing_count, modified_tree) if all gene copies are unique:
        missing_count (Integer): Number of taxa missing from the tree.
        modified_tree (String): Modified Newick string with simplified taxa labels.
    -Returns Nothing if there are repeated gene copies.

Examples: 
# Case 1: All taxa appears for once but there are F, G, H missing (number of missing = 3). 
Input: count_missing("(A_0_1:0.01,(((B_0_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
Output: (3, "(A_1:0.01,(((B_1:0.004,C_1:0.003):0.007,(D_0:0.02,E_1:0.02):0.01)));")

# Case 2: # B, G, H missing (number of missing = 3) and there are repeated copies A_0_1 and A_1_1. Returns nothing
Input: count_missing("(F_0_1:0.8,(A_0_1:0.01,((A_1_1:0.004,C_1_1:0.003):0.007,(D_1_0:0.02,E_0_1:0.02):0.01)));")
Output: Nothing 
""" 
function count_missing(newick::String)
    taxa = "ABCDEFGH"  # Species tree only have those taxa
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
Input: Input files with newick trees from simphy
Output: Modified strings got output to an output files 
- If the input newick has repeated taxa, not included in the downstream analysis
- If the input newick has more than 4 taxa missing (three taxa left), not included in the downstream analysis
Notes: If the input has repeated taxa, then paralohy is not hidden. To simulate hidden paralogy, then the tree will not be processed. 
If no hidden paralogy got simulated, then dup_rate = 0 and loss_rate = 0, so there is no gene trees got removed. However, this code will still re-name the tree tips so that A_locusID_accessionID is changed to A_accessionID. 
"""
function modify_newicks(input::String, output::String)
    lines = readlines(input)
    valid_lines = String[]  
    
    for nwk in lines
        temp = count_missing(nwk)
        if isnothing(temp)
            println("The input $input has repeated taxa, not proceeded to downstream analysis")
            continue
        elseif temp[1] > 4
            println("The input $input has more than four taxa missing, not proceeded to downstream analysis")
            continue
        else
            push!(valid_lines, temp[2] * "\n")
            # println("Successfully modify $input !")
        end
    end
    
    if !isempty(valid_lines) # Only write the file if there's valid data to output
        open(output, "w") do io
            write(io, join(valid_lines))
        end
    else
        println("No valid Newick strings were found for $input. No output file created.")
    end
end

"""
Run modify_newicks given the number of gene trees (newicks) 
Input: 
    -input: input dir path 
    -output: output dir path 
    -n: Number of gene trees / newicks to be modified
"""
function modify_newick_for_n_genes(input::String, output::String, n::Int, iteration::Int)
    for gene_tree in 1:n
        genenum_string = pad_number(gene_tree, n) # see Part 3 "pad_number" below
        output_tree_name = string("g_trees_noLocusID_Gene$(genenum_string)_Int$(iteration).trees")
        input_dir = joinpath(input, "g_trees$genenum_string.trees")
        output_dir = joinpath(output, output_tree_name)
        modify_newicks(input_dir, output_dir)
    end
end 

#-----------------------------------------------#       
#               Part 2 
#    Re-run simphy until:
# 1) n_genes hit the target; or 
# 2) max_iteration reaches   
#-----------------------------------------------#
""" 
Writes N random integers to be used as seeds in imphy simulation
Input:  -master_seed: An Int to be used in Xoshiro() to generate random num 
        -n_reps: number of reps to be simulated
        -max_iteration: maximum iteration to simulate simphy
        -output_dir: The output directory to store the random_seeds.txt file 
Output: - A random_seed.txt file with n_reps x max_iterations random seeds in the output_dir
            First line: master seed
            Second line: list all iterations
            Lines below: list n_rep and the seeds 
        - Return an array with all seeds
""" 
function seed_generator(master_seed::Int, n::Int, m::Int, output_dir::String, output_file_name::String) 

    rng = StableRNG(master_seed)
    seeds = abs.(rand(rng, Int32, n, m)) # generates n x m seed array with positive numbers only 

    # write all seeds to a random_seeds.txt file in the output_dir 
    output_file = joinpath(output_dir, output_file_name)
    rep = 1
    open(output_file, "w") do io
        for j in 1:m
            print(io, j == m ? "col$j\n" : "col$j\t") # if j is the last then separate by a new line
        end

        for i in 1:n
            for j in 1:m
                print(io, j == m ? "$(seeds[i, j])\n" : "$(seeds[i, j])\t") # if j is the last then separate by a new line
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
    -int: Rounded-up estimated batch size for the next simulation. 
Calculation explanations:
    For a missing rate of 0.2 and 30 target genes:
    - Missing trees: 30 * 0.2 = 6
    - Success rate: 1 - 0.2 = 0.8
    - estimated batch size to be simulated: 6 / 0.8 = 7.5 → 8 (rounded up)
    - If success rate is 0 (cannot be divided), so estimated batch = n_genes 
"""
function calculate_batch(num::Int, n_genes::Int) 
    num_missing = n_genes - num 
    success_rate = num / n_genes 
    if success_rate == 0
        estimated_batch = n_genes 
    else 
        estimated_batch = round(num_missing / success_rate) 
    end 
    # println("For replicate", n_reps, " estimated batch size = ", estimated_batch)
    return Int(estimated_batch)
end 

"""
pad_number(num::Int, range::Int) -> String
Converts `num` into a zero-padded string based on the number of digits required 
by `range`. This function is required because simphy outputs number as below. 
Inouts: 
    - `num::Int`: The number to pad.
    - `range::Int`: The range determining the padding width.
Outputs: 
    - `String`: A zero-padded string representation of `num`.
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

#-----------------------------------------------#       
#               Part 3  
#       Other utilities functions     
#-----------------------------------------------#
"""
Checks if a specified path (`dir`) exists and rm or keep the files based on user inputs 
    The function is useful for testing the same parameters repeatedly. 
- If it exists, the function allows the user to decide whether to remove it. 
- If there are more than one paths, user input "ALL" (case-insensitive) to remove all paths 
- If there are one path, ->`y` (case-insensitive), the directory or file is removed using `rm`.
- If the user inputs `n` (case-insensitive), the function raises an error and stops further execution.
- If the user inputs anything else, the function repeatedly prompts for valid input (`y` or `n`).
- If the path does not exist, the function is skipped 
"""
function check_existing_dir(dir_input::Vector)

    len = length(dir_input)
    if len == 0 # If input an empty string, then output an error 
        error("Oh no! It's empty string/vector. Check out the code.")
        return 
    end

    if len > 1 && ispath(dir_input[1]) # if input is a list of strings and the first path of the list exists 
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
This function replaces text tips with letters (max # taxa: 26). Examples:

1. String input → string output:  
Input: replace_tips_with_letters("(Homo:0.01,(((Cat:0.004,Dog:0.003):0.007,(Dinosaur:0.02,Birds:0.02):0.01)));")  
Output: "(A:0.01,(((B:0.004,C:0.003):0.007,(D:0.02,E:0.02):0.01)));"

2. Tree input → tree output:  
Input: replace_tips_with_letters(readTopology("(Homo:0.01,(((Cat:0.004,Dog:0.003):0.007,(Dinosaur:0.02,Birds:0.02):0.01)));"))  
Output:  
HybridNetwork, Rooted Network  
9 edges, 10 nodes (5 tips, 5 internal).  
Tip labels: A, B, C, D, ...  
"(A:0.01,(((B:0.004,C:0.003):0.007,(D:0.02,E:0.02):0.01)));"

3. Handles special characters like *:  
Input: replace_tips_with_letters("(Homo:0.01*4,(((Cat:0.004*5,Dog:0.003*10):0.007,(Dinosaur:0.02*20,Birds:0.02*3):0.01*9)));")  
Output: "(A:0.01*4,(((B:0.004*5,C:0.003*10):0.007,(D:0.02*20,E:0.02*3):0.01*9)));"
"""
function replace_tips_with_letters(tree::Union{HybridNetwork, String}) 
	#= Goal: replace species tips with A, B, C, etc, but tips should be less than or equal to 26. 
    The input (tree) can both be a string or a tree
    If input is tree, then return a tree object
    If input is a string, then return a string: which could be used to modified the input to SimPhy=#
    new_labels = ['A':'Z';]
    if isa(tree, HybridNetwork) # If input is a tree, output a tree 
        tips = tipLabels(tree)

        if length(tips) > 26 # Give an warning and return the original tree if there are more than 26 tips 
            println("Warning: The number of tips exceeds 26. That's too much. Exiting.")
            return tree
        end 
    
        for (i, tip) in enumerate(tree.leaf)
            tip.name = string(new_labels[i])
        end
        return tree
    
    elseif isa(tree, String) # if input is a Newick string, output a string 
        tips = split(tree, [',', '(', ')', ':', '*', ';'])  # Extract tips
        tips = filter(x -> !(x in ["", ":", ";"]) && !occursin(r"\d", x), tips) # Remove empty, :, and any string containing numbers

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
Return output path string based on inputs: 
  -folder_path_list: A list containing folder path from rep 1 to n_reps
  -simulation_rep: Current rep ID
  -path_string: The appending string to extend the folder name
"""
function setup_rep_output_folders(folder_path_list::Vector, simulation_rep::Int, path_string::String)
  rep_folder_path = folder_path_list[simulation_rep]  
  return (joinpath(rep_folder_path, path_string))
end 

"""
A test function which is not used in the main script
    Input: 
        -file_path: A file path to a tree file containing multiple boostrap trees 
        -Output_dir: The output directory to store resulting files as below 
        -output_filename: a string to 
        -relative_path: Boolean -> Create a relative path or absolte in the output filename 
            -> If true -> create relative path else absolute path 
    Return: 
        -Split the each individual tree in the file path into a separate file
        -Write the file path to each individual tree into a .txt file 
This function is not used in the major script but used as a test function to test if the input for bootsnaq could be converted as [[Hybird tree]]
"""
function write_boostrapTrees_2_Ind_filepath(file_path::String, output_dir::String, output_filename::String, relative_path::Bool)
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
A function to replace occurrences of "POP" in a file with the corresponding letter extracted from the first column.
    Input: 
        - input_file: File path to a .ind file where the third column may contain "POP"
        - output_file: The output file path
    Return: 
        - A new file with updated values in the third column 
            EX: If first col "A_0", then the assigned pop is A 
                If first col "A_1", then the assigned pop is also A 
            This is because A_0 and A_1 are accessions from the same pop
This function is used to process .ind input files where "POP" is a placeholder for group labels based on the first column.
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
Generate a dictionary of stable random seeds for each software name using a given master seed.
Inputs:
    - master_seed(Int): The initial master seed
    - software_names(Vector{String}): A list of software names that require unique seeds.
Returns
    A dictionary where keys are software names and values are stable integer seeds.
Example: 
    master_seed = 123456
    software_names = ["SimPhy", "RAxML", "Astral", "SeqGen"]
    software_seeds = generate_software_seeds(master_seed, software_names)
"""
function generate_software_seeds(master_seed::Int, software_names::Vector{String})
    rng = StableRNG(master_seed) 
    seeds = abs.(rand(rng, Int, length(software_names)))  
    return Dict(software_names[i] => seeds[i] for i in eachindex(software_names))
end

"""
generate_master_seed(params::Dict{String, Int})
Generate a master seed based on given parameter values. The parameters can be:
    - continuous values (e.g., 0.01, 0.0001) → encoded by assigning digits after the decimal point to primes.
    - categorical values → mapped manually using `category_map` (hard coded inside the function)
Returns an `Int` seed, which is the product of prime values corresponding to the parameter settings.
"""
function generate_master_seed(params::Dict{String, Any})

    digit_to_prime = Dict('0'=>2, '1'=>3, '2'=>5, '3'=>7, '4'=> 11, '5'=>13, '6'=>17, '7'=>19, '8'=>23, '9'=>29)
    category_map = Dict("ratevar" => Dict("N"=>37, "G"=>41, "LG"=>43, "GL" => 43, "G*L" => 43))
    master_seed = big(1)  # big to prevent large integer overflow 
    prime = 1 # prime starts from 1 

    for (param, value) in params
        
        if haskey(category_map, param) # If the parameter is in category_map 
            prime = category_map[param][value]

        elseif isa(value, Float64) || isa(value, Float32) # If the parameter takes a number   
            splitted_number = split(string(value), '.') # split the string based on the decimal dot

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
            prime = nextprime(value) # nextprime assigns a prime number which is bigger than the value
        else
            error("Oh no something is wrong with the parameters setting. Unable to generate master seed")
        end   
        master_seed *= prime
    end
    return Int(master_seed) # change to integer because generate_software_seeds only takes integers 
end

"""
Function to get the dictionary for setting up the master seed
    Params: The parameters are based on the string "paramname_root", see simulation_iqtree.jl 
        Note: The parts and params are hard-coded. 
    Return: output will be directly used as the input for the function generate_master_seed(params)
"""
function get_dict_for_seed_setting(params_string::String) 

    parts = split(params_string, "-")

    dup_rate = parse(Float32, split(parts[1], "DUP")[2])
    loss_rate = parse(Float32, split(parts[1], "DUP")[2])
    ratevar = split(parts[3],"RV")[2]
    N_ind = parse(Int, split(parts[4],"N_ind")[2])

    params = Dict{String, Any}(
           "dup_rate" => dup_rate,
           "LOSS_RATE" => loss_rate,
           "ratevar" => ratevar, 
           "n_inds" => N_ind)

    return params 
end 