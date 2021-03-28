function BestHit(leaves::Dict{Taxon,Int})
    return last(findmax(leaves))
end

function LCA(leaves::Dict{Taxon,BlastResult},minimal::Float64)
    besthitscore = first(findmax([last(l).bitscore for l in leaves]))
    filter!(x -> last(x).bitscore < besthitscore*minimal, leaves)
    taxa = collect(keys(leaves))
    return lca_node = lca(taxa)
end

function weightedLCA(leaves::Dict{Taxon,Int}, minimal::Float64, cutoff::Float64)
    @assert cutoff > 0.5 && cutoff < 1
    besthitscore = first(findmax([last(l).bitscore for l in leaves]))
    filter!(x -> last(x).bitscore < besthitscore*minimal, leaves)

    taxa = collect(keys(leaves))
    tree = topolgoy(taxa)

    total_bitscore = sum([last(l).bitscore for l in leaves])
    threshold_bitscore = cutoff * total_bitscore

    next = Stack{PhyloTree}()
    push!(next,tree)
    current_lca = tree
    while true
        current_node = pop!(next)
        sub_bitscore = sum([leaves[leave.node].bitscore for leave in Leaves(tree)])
        if sub_bitscore >= threshold_bitscore
            current_lca = current_node
            children_tree =  current_node.children
            for ctree in children_tree; push!(next,ctree); end
        end
        if isempty(next)
            break
        end
    end
    return current_lca.node
end