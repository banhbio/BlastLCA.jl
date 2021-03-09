function blastLCA(;filepath::AbstractString, outpath::AbstractString, sqlite::SQLite.DB, taxonomy::Taxonomy.DB, method::Function, header::Bool=false)
    f = open(filepath,"r")
    o = open(outpath,"w")
    blastLCA(f, o; sqlite=sqlite, taxonomy=taxonomy, method=method, header=header)
    close(f)
    close(o)
end

function blastLCA(f::IOStream, o::IOStream; sqlite::SQLite.DB, taxonomy::Taxonomy.DB, method::Function, header::Bool=false)
    header ? readline(f) : nothing #skip header

    bitscores = Dict{Taxon,Int}()

    line = readline(f) #read first line
    while true
        rows = split(line, "\t")    
        @assert length(rows) == 12
        qseqid = rows[1]
        sseqid = rows[2]
        bitscore = rows[12]
        taxid = get(sqlite, sseqid, nothing)
        @assert taxid !== nothing

        taxon = Taxon(taxid,taxonomy)
        if ! haskey(bitscores, taxon)
            bitscores[taxon] = bitscore
        else
            if bitscore > bitscores[taxon]
                bitscores[taxon] = bitscore
            end
        end

        next = readline(f)
        next_rows = split(line, "\t")
        next_qseqid = next_rows[1]

        if qseqid != next_qseqid
            assignment = method(bitscores)
            write(o,"$(qseqid)\t$(assignment)")
            bitscores = Dict{Taxon,Int}() #initialize
        end

        line = next
        
        eof(f) ? break : nothing
    end
end