module BlastLCA

using Taxonomy
using SQLite
export SQLite
using DataStructures:
    Stack
export create!,insert!,get,
       BlastResult,
       blastLCA,
       BestHit,freeLCA,weightedLCA
include("main.jl")
include("lca.jl")
include("database.jl")
end