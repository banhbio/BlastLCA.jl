function freeLCA(input::Tuple{String,Dict{Taxon,BlastResult}}, minimal::Float64,ranks::Vector{Symbol},precision::Dict{Symbol, Float64})
    leaves = last(input)
    besthitscore = first(findmax([last(l).bitscore for l in leaves]))
    filter!(x -> last(x).bitscore > besthitscore*minimal, leaves)
    taxa = collect(keys(leaves))
    current_lca = topolgoy(taxa)
    corrected_taxon = cut_by_precision(current_lca, ranks, precision, leaves)
    return corrected_taxon
end

function weightedLCA(input::Tuple{String,Dict{Taxon,BlastResult}}, minimal::Float64, cutoff::Float64, ranks::Vector{Symbol},precision::Dict{Symbol, Float64})
    @assert cutoff > 0.5 && cutoff < 1
    leaves = last(input)
    besthitscore = first(findmax([last(l).bitscore for l in leaves]))
    filtered_leaves = filter(x -> last(x).bitscore > besthitscore*minimal, leaves)

    taxa = collect(keys(filtered_leaves))
    tree = topolgoy(taxa)

    total_bitscore = sum([last(l).bitscore for l in filtered_leaves])
    threshold_bitscore = cutoff * total_bitscore

    next = Stack{PhyloTree}()
    push!(next,tree)
    current_lca = tree
    while true
        current_node = pop!(next)
        sub_bitscore = sum([leaves[leave.node].bitscore for leave in Leaves(current_node)])
        if sub_bitscore >= threshold_bitscore
            current_lca = current_node
            children_tree =  current_node.children
            for ctree in children_tree; push!(next,ctree); end
        end
        if isempty(next)
            break
        end
    end
    corrected_taxon = cut_by_precision(current_lca, ranks, precision, leaves)
    return corrected_taxon
end

function cut_by_precision(current_lca::PhyloTree, ranks::Vector{Symbol}, precision::Dict{Symbol, Float64}, leaves::Dict{Taxon,BlastResult})
    max_sub_pident = maximum([leaves[leave.node].pident for leave in Leaves(current_lca)])
    lineage = Lineage(current_lca.node)
    for taxon in lineage
        r = rank(taxon)
        if ! (r in ranks) || !haskey(precision, r)
            continue
        end
        if precision[r] > max_sub_pident
            return taxon
        end
    end
    return lineage[end]
end