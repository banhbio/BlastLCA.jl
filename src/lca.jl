function BestHit(leaves::Dict{Taxon,Int})
    return last(findmax(leaves))
end

function freeLCA(leaves::Dict{Taxon,BlastResult},minimal::Float64,ranks::Vector{Symbol},precision::Dict{Symbol, Float64})
    besthitscore = first(findmax([last(l).bitscore for l in leaves]))
    filter!(x -> last(x).bitscore < besthitscore*minimal, leaves)
    taxa = collect(keys(leaves))
    taxon = lca(taxa)
    lineage = Lineage(taxon)
    reformated_lineage = cut_by_precision(lineage, ranks, precision)
    return reformated_lineage
end

function weightedLCA(leaves::Dict{Taxon,BlastResult}, minimal::Float64, cutoff::Float64, ranks::Vector{Symbol},precision::Dict{Symbol, Float64})
    @assert cutoff > 0.5 && cutoff < 1
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
    lineage = Lineage(current_lca.node)
    corrected_rank = cut_by_precision(current_lca, ranks, precision, leaves)
    corrected_lineage = reformat(lineage, ranks)[Until(corrected_rank)]
    return corrected_lineage
end

function cut_by_precision(current_lca::PhyloTree, ranks::Vector{Symbol}, precision::Dict{Symbol, Float64}, leaves::Dict{Taxon,BlastResult})
    max_sub_pident = maximum([leaves[leave.node].pident for leave in Leaves(current_lca)])
    for r in ranks
        if !haskey(precision, r)
            continue
        end
        if precision[r] > max_sub_pident
            return r
        end
    end
    return ranks[end]
end