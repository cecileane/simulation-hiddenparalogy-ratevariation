using PhyloNetworks

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






