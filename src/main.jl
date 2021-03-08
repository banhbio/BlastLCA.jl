function blastLCA(;filepath::AbstractString, outpath::AbstractString, sqlite::SQLite.DB, taxonomy::Taxonomy.DB, method::Function, header::Bool=false)
    f=open(filepath,"r")
    g=open(outpath,"w")

    header ? readline(f) : nothing

    current_qseqid = ""
    bitscores = Dict{Taxon,Int}()

    while !eof(f)
        line = readline(f)
        rows = split(line, "\t")
        
        @assert length(row) == 12

        qseqid = rows[1]

        sseqid = rows[2]
        bitscore = rows[12]

        taxid = get(sqlite, sseqid, nothing)
        @assert taxid !== nothing

        taxon = Taxon(taxid,taxonomydb)
        if ! haskey(d, taxon)
            bitscores[taxon] = bitscore
        else
            if bitscore > d[taxon]
                bitscores[taxon] = bitscore
            end
        end

        next = readline(f)
        next_rows = split(line, "\t")
        next_qseqid = next_rows[1]

        if qseqid != next_qseqid
            assignment = method(bitscores)
            write(g,"$(qseqid)\t$(assignment)")
            qseqid = next_qseqid
            bitscores = Dict{Taxon,Int}()
        end

        rows = next_rows
    end
    close(f)
    close(g)
end