module BlastLCA

using Taxonomy
using SQLite
using CSV
using AbstractTrees
using DataStructures:
    Stack
using DataFrames

export BlastResult,
       blastLCA,
       BestHit,freeLCA,weightedLCA
include("main.jl")
include("tree.jl")
include("lca.jl")
end