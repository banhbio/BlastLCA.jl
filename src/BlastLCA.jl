module BlastLCA

using Taxonomy
using SQLite
export SQLite

include("main.jl")
include("lca.jl")
include("database.jl")
include("struct.jl")

end
