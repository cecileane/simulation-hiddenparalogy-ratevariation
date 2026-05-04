# The script to run snaq 
using Distributed
using ArgParse
@everywhere using CSV
@everywhere using DataFrames
@everywhere using QuartetNetworkGoodnessFit
@everywhere using PhyloNetworks
@everywhere using SNaQ
@everywhere include("utilities.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        # Folders passed down from simulation.jl for meaningful snaq filenames

        "--outfolder" 
        help = "root folder to be specified"
        arg_type = String
        required = true

        "--iqtreefolder"
        help = "iqtree output folder"
        arg_type = String
        required = true

        "--astralfolder"
        help = "astral output folder"
        arg_type = String
        required = true

        "--snaqfolder"
        help = "snaq output folder"
        arg_type = String
        required = true

        "--paramname_root" # create meaningful file names for snaq outputs
        help = "Parameter set used to simulate sim-phy"
        arg_type = String
        required = true

        "--seed_net0"
        help = "Seed for estimating tree Hmax = 0"
        arg_type = Int
        required = true

        "--seed_net1"
        help = "Seed for estimating network Hmax = 1"
        arg_type = Int
        required = true 

        "--seed_qgof0"
        help = "Seed for QuartetNetworkGoodnessFit for Hmax = 0"
        arg_type = Int
        required = true

        "--seed_qgof1"
        help = "Seed for QuartetNetworkGoodnessFit for Hmax = 1"
        arg_type = Int
        required = true

        "--seed_net0_bs"
        help = "bootnet seed Hmax = 0"
        arg_type = Int 
        required = true 

        "--seed_net1_bs"
        help = "bootnet seed Hmax = 1" 
        arg_type = Int
        required = true

        "--runs"
        help = "Number of runs to run SnaQ"
        arg_type = Int
        required = true

        "--n_snaqboot_rep"
        help = "Number of replicate for SnaQ boostrap for Hmax = 1" 
        arg_type = Int
        required = true

        "--n_inds"
        help = "Number of individuals per species" 
        arg_type = Int
        required = true

        "--method"
        help = """Method to run sanq pipeline. Options are:
                1. 'snaq_bootstrap' to run boostrapping of snaq
                2. 'QuartetNetworkGoodnessFit' to compare goodness of fit"""
        arg_type = String
        required = true # default is set in snaq.jl 

    end 
    return parse_args(s)
end

#-------------------------------------# 
# Parse commandline arguments 
#-------------------------------------# 
parsed_args = parse_commandline()

# Get arguments: 
outfolder = parsed_args["outfolder"]
iqtreefolder = parsed_args["iqtreefolder"]
astralfolder = parsed_args["astralfolder"]
snaqfolder = parsed_args["snaqfolder"]
paramname_root = parsed_args["paramname_root"]

# Get args used for SnaQ: 
runs = parsed_args["runs"]
n_snaqboot_rep = parsed_args["n_snaqboot_rep"]
seed_net0 = parsed_args["seed_net0"]
seed_net0_bs = parsed_args["seed_net0_bs"]
seed_net1 = parsed_args["seed_net1"]
seed_net1_bs = parsed_args["seed_net1_bs"]
n_inds = parsed_args["n_inds"]
seed_qgof0 = parsed_args["seed_qgof0"]
seed_qgof1 = parsed_args["seed_qgof1"] 

# The methods to run SNAQ: 
# TWO options: 
# 1. run boostrapping of snaq 
# 2. run snaq and then compare loglikelihoods 
method = parsed_args["method"] 
# method = "QuartetNetworkGoodnessFit" # or "snaq_bootstrap"
if !(method in ["snaq_bootstrap", "quartetnetworkgoodnessfit", "both"])
    error("method must be 'snaq_bootstrap', 'quartetnetworkgoodnessfit'," *
        " or 'both'")
end

# Specify folder path: 
species_tree = joinpath(astralfolder,"astral.tre") # Astral tree 
gene_trees = joinpath(iqtreefolder,"besttrees.tre") # IQtree trees
outdir = snaqfolder # output directory for snaq. Each rep gets its own folder

# Make output directories for H=0 and H=1
H0folder = joinpath(outdir, "H0_output")
H1folder = joinpath(outdir, "H1_output")
mkpath(H0folder)
mkpath(H1folder)

# Load starting topology, gene trees, and CFs
# species_tree = readMultiTopology(species_tree)[102] 
# Older astral (5.7.8): starting topology is the last (102nd) tree
# Newer astral: starting topology is the first tree
# Read the astral tree file and clean NaN values before parsing
species_tree_content = read(species_tree, String)
species_tree_content = clean_newick_nan(species_tree_content)
species_tree = readnewick(species_tree_content)
# readnewick from PhyloNetworks 1.1.0 

# Now, the tree tip is "A", "B", "C", etc, instead of "A_0", "A_1", etc. 
gene_trees = readmultinewick(gene_trees)
taxonmap, mappingfile = map_accessions_to_species_dict(gene_trees, outdir) 
# The above function retuns a Dict with species => individuals 
# id_to_species.csv saved in outdir (snaqfolder), unique per replicate

#-------------------------------------#
# Calculate concordance factors
#-------------------------------------# 

if n_inds == 1 # Probably multiple gene copies per species

    # when n_inds == 1, we don't need to have taxonmap
    # but this is just for consistency and changing tipnames to A, B, C, etc 
    # see https://juliaphylo.github.io/SNaQ.jl/dev/man/snaq_est/

    q,t = countquartetsintrees(gene_trees, taxonmap) 
    nt = tablequartetCF(q,t) 
    df = DataFrame(nt, copycols=false) 
    CSV.write("$outdir/CF_results.csv", df) # save CFs into a file 
    iqtreeCF = readtableCF(df)

    num_cf = nrow(df) 

else # multiple individuals per species  
    q,t = countquartetsintrees(gene_trees)
    df_ind = DataFrame(tablequartetCF(q,t)) 
    # df_ind is the table of CF at individual level 
    ind_cf_path = joinpath(outdir, "tableCF_individuals.csv")  
    CSV.write(ind_cf_path, df_ind) # save CFs into a file
    df_sp = mapallelesCFtable(mappingfile, 
                            ind_cf_path;
                            columns=2:5)
    d_sp= readtableCF!(df_sp, mergerows=true) # DataCF object 

    # Omit external branch length estimation (slow)
    # Check the documenation for this section: 
    # https://juliaphylo.github.io/SNaQ.jl/dev/man/multiplealleles/ 
    # -> calculated by averaging the CFs of quartets of individuals 
    df_sp_ave = DataFrame(tablequartetCF(d_sp)) 
    df_sp_reduced = filter(!hasrep, df_sp_ave)  
    # Removes repeated-taxa rows; terminal branch lengths not estimated
    CSV.write("$outdir/CF_results.csv", df_sp_reduced)
    iqtreeCF = readtableCF(df_sp_reduced) 

    num_cf = nrow(df_sp_reduced) 
    # num_cf used in capushe model comparison as number of data points
end 

#-------------------------------------#
# Run SNaQ
#-------------------------------------# 
#= 
Run SNaQ: 
Option 1: boostrapping of snaq
    -> run snaq with Hmax = 0 and Hmax = 1
    -> run bootsnaq with Hmax = 0 and Hmax = 1
    The confidence of hybird edge might show useful information 
    However, this cannot help us compare between networks with different Hmax
    because boostrapping assumes that Hmax is fixed 

Option 2: run snaq and then compare loglikelihoods  
    -> run snaq with Hmax = 0, Hmax = 1, Hmax = 2, etc.
    -> compare the loglikelihoods of the two networks 
=# 

# addprocs(processors)
# hmax = 0 --> The starting tree for Hmax = 1 
# Below estimates network 
# Here snaq for H = 0 does not need that many runs to find the optimal tree
# because the starting tree is already a good estimate from astral
# Using the default 10 runs is fine 
println("SNaQ Hmax=0: runs=10, seed=$seed_net0. Running...")
net0out = joinpath(H0folder, "H0")
# Testing: runs_net0 = runs; production (runs > 50): runs_net0 = 10
runs_net0 = runs > 50 ? 10 : runs
net0 = snaq!(species_tree, iqtreeCF, hmax=0, filename=net0out,
    seed=seed_net0, runs=runs_net0)

# hmax = 1 --> The network with one hybrid edge 
println("SNaQ Hmax=1: runs=$runs, seed=$seed_net1. Running...")
net1out = joinpath(H1folder, "H1")
net1 = snaq!(net0, iqtreeCF, hmax=1, filename=net1out,
    seed=seed_net1, runs=runs)

#-------------------------------------#
# Post-SNaQ analyses 
# --> Extract gamma values
# --> Check if we get the true species tree for Hmax = 0 
# --> Check if we get the true species tree displayed for Hmax = 1 
# --> boostrapping of snaq (optional)
# --> Compare the goodness of fit between net0 and net1
#-------------------------------------#  

# Extract gamma values, saved in to summary csv file 
# This will only be applied to H1
gammas = [e.gamma for e in net1.edge if e.hybrid]  
println("Number of hybrid edges found: ", length(gammas))
println("Gamma values: ", gammas)

# With hmax=1, we expect exactly 2 gamma values. 
if length(gammas) >= 2
    gamma_1 = gammas[1] 
    gamma_2 = gammas[2] 
    if gamma_1 < gamma_2 # Ensure gamma_1 is always bigger than gamma_2 
        gamma_1, gamma_2 = gamma_2, gamma_1
    end
else
    error("Expected ≥2 hybrid edges with hmax=1, found $(length(gammas))")
end 

# Run boostrapping only if specified 
if method == "snaq_bootstrap" || method == "both"
    # Legacy code for boostrapping
    # We will most likly not use this

    # Bootstrapping
    println("bootsnaq Hmax=0: runs=$runs, seed=$seed_net0_bs. Running...")
    bslist = joinpath(iqtreefolder, "bslist.txt")
    gene_tree_list = [ [tree] for tree in gene_trees ]

    for tree_wrapper in gene_tree_list # simplify tip labels to species names
        tree = tree_wrapper[1]
        simplified_tree = simplify_tip_labels(tree)
        tree_wrapper[1] = simplified_tree
    end # One tip per species

    # Boostrapping H = 0:  
    bootnet0 = bootsnaq(
        species_tree,
        gene_tree_list,
        hmax=0,
        nrep=n_snaqboot_rep,
        filename=net1out,
        seed=seed_net0_bs,
        runs=runs)
    net0 = readnewick("$net0out.out")
    BSe_tree0, tree0 = treeEdgesBootstrap(bootnet0,net0)
    bootnet0_output = joinpath(H0folder,"bootnet0_onTree0.csv")
    CSV.write(bootnet0_output, BSe_tree0)

    # Boostrapping H = 1: 
    println("bootsnaq Hmax=1: runs=$runs, seed=$seed_net1_bs. Running...")
    bootnet1 = bootsnaq(
        net0,
        gene_tree_list,
        hmax=1,
        nrep=n_snaqboot_rep,
        filename=net1out,
        seed=seed_net1_bs,
        runs=runs)
    net1 = readnewick("$net1out.out")
    BSe_tree1, tree1 = treeEdgesBootstrap(bootnet1,net1)
    bootnet1_output = joinpath(H1folder,"bootnet1_onNet1.csv")

    CSV.write(bootnet1_output, BSe_tree1)

end 

if method == "quartetnetworkgoodnessfit" || method == "both" 

    # Compare the goodness of fit between net0 and net1 
    # see tutorial: 
    # see QuartetNetworkGoodnessFit.jl docs: /stable/man/gof/#goodness_of_fit_1
    
    # load CFs from the saved csv file 
    cffile = joinpath(outdir, "CF_results.csv") 
    qCF = CSV.read(cffile, DataFrame)

    # Check for quartets with zero gene trees
    # quarnetGoFtest! cannot handle missing data (ngenes = 0)
    zero_gene_quartets = filter(row -> row.ngenes == 0, qCF)
    
    if nrow(zero_gene_quartets) > 0
        # Skip goodness-of-fit test if there are missing quartets
        println("\n" * "="^70)
        println("WARNING: $(nrow(zero_gene_quartets)) quartet(s) with ngenes=0")
        println("This occurs due to stochastic gene loss in the simulation.")
        println("Skipping quarnetGoFtest for this replicate.")
        println("Missing quartets:")
        for row in eachrow(zero_gene_quartets)
            println("  $(row.t1), $(row.t2), $(row.t3), $(row.t4)")
        end
        println("="^70 * "\n")
        
        # Still save network information and gamma values without GoF test
        score_net0 = loglik(net0)
        score_net1 = loglik(net1)
        
        # Create minimal output files indicating the test was skipped
        gof0_path = joinpath(H0folder, "snaq_gof_results_H0.csv")
        gof1_path = joinpath(H1folder, "snaq_gof_results_H1.csv")
        
        df_0 = DataFrame(type = ["score", "status"], 
                        value = [score_net0, "SKIPPED: missing quartets"])
        CSV.write(gof0_path, df_0)
        
        df_1 = DataFrame(type = ["score", "gamma_1", "gamma_2", "status"], 
                        value = [score_net1, gamma_1, gamma_2,
                            "SKIPPED: missing quartets"])
        CSV.write(gof1_path, df_1)
        
        println("Saved minimal results (without GoF test) to:")
        println("  $gof0_path")
        println("  $gof1_path")
        
    else
        # Proceed with normal goodness-of-fit test
        score_net0 = loglik(net0)
        score_net1 = loglik(net1)

        # This pipeline is based on: 
        # see QuartetNetworkGoodnessFit.jl docs: /stable/man/gof/
        # Warm-up with few sims to find optimal branch lengths, then full run
        res0 = quarnetGoFtest!(net0, qCF, true; seed=201, nsim=10); #warm-up 
        res1 = quarnetGoFtest!(net1, qCF, true; seed=202, nsim=10); #warm-up 
        net0 = res0[5] 
        net1 = res1[5] 
        res0 = quarnetGoFtest!(net0, qCF, false; seed=seed_qgof0, nsim=1000);
        res1 = quarnetGoFtest!(net1, qCF, false; seed=seed_qgof1, nsim=1000);
        # res0 = quarnetGoFtest!(net0, qCF, true; seed=seed_qgof0, nsim=1000);
        # res1 = quarnetGoFtest!(net1, qCF, true; seed=seed_qgof1, nsim=1000);
        # Branch length difference between true/false is small but noted

        # Save the results
        gof0_path = joinpath(H0folder, "snaq_gof_results_H0.csv")
        gof1_path = joinpath(H1folder, "snaq_gof_results_H1.csv")

        # save the results for H = 0 
        # Define meaningful names for each result element from quarnetGoFtest!
        # Names based on QuartetNetworkGoodnessFit.jl docs and observed output
        result_names = ["p", "z_uncorrected", "sigma", "bootstrap_values", 
                        "network", "z_bootstrap"]
        
        serializable_res0 = []
        for (i, x) in enumerate(res0)
            name = i <= length(result_names) ? result_names[i] : "unknown_$i"
            if typeof(x) <: Number
                push!(serializable_res0, (type=name, value=x))
            elseif typeof(x) <: AbstractArray
                push!(serializable_res0, (type=name, value=join(x, ";")))
            elseif typeof(x) == HybridNetwork
                push!(serializable_res0, (type=name, value=string(x)))
            else
                push!(serializable_res0, (type=name, value=string(x)))
            end
        end
        
        # build net0 summary dataframe 
        df_0 = DataFrame(serializable_res0)
        # Add loglik as a separate row for H = 0
        loglik_row = DataFrame(type = ["score"], value = [score_net0])
        df_0 = vcat(df_0, loglik_row)

        println(df_0) # print to check 

        CSV.write(gof0_path, df_0)

        # save the results for H = 1 
        serializable_res1 = []
        for (i, x) in enumerate(res1)
            name = i <= length(result_names) ? result_names[i] : "unknown_$i"
            if typeof(x) <: Number
                push!(serializable_res1, (type=name, value=x))
            elseif typeof(x) <: AbstractArray
                push!(serializable_res1, (type=name, value=join(x, ";")))
            elseif typeof(x) == HybridNetwork
                push!(serializable_res1, (type=name, value=string(x)))
            else
                push!(serializable_res1, (type=name, value=string(x)))
            end
        end

        # build net1 summary dataframe 
        df_1 = DataFrame(serializable_res1)
        # Add loglik as a separate row for H = 1
        loglik_row = DataFrame(type = ["score"], value = [score_net1])
        df_1 = vcat(df_1, loglik_row)

        # save gamma values to gof1 output
        gamma_rows = DataFrame(type = ["gamma_1", "gamma_2"], 
                            value = [gamma_1, gamma_2])
        df_1 = vcat(df_1, gamma_rows)

        println(df_1) # print to check

        CSV.write(gof1_path, df_1)
    
    end  # end of if/else for missing quartets check 

end





