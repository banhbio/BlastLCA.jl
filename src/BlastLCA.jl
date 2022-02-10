module BlastLCA

using Taxonomy
using SQLite
using CSV
using DataStructures:
    Stack
using DataFrames
export create!,insert!,get,
       BlastResult,
       blastLCA,
       BestHit,freeLCA,weightedLCA
include("main.jl")
include("lca.jl")
include("database.jl")
end