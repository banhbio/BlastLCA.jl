function blastLCA(;filepath::AbstractString, outpath::AbstractString, blastLCAdb::BlastLCA.DB, taxonomydb::Taxonomy.DB, evalue::Float64, minimal::Float64,cutoff::Float64,header::Bool=false)
    f=open(filepath,"r")
    g=open(outpath,"w")

    header ? readline(f): nothing

    current_qseqid = ""
    bitscores = Dict{Taxon,Int}()
    line = readline(f)
    rows = split(line, "\t")
    while true
        @assert length(row) == 12

        qseqid = rows[1]

        sseqid = rows[2]
        bitscore = rows[12]

        taxid = get(blastLCAdb, sseqid, nothing)
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
            assignment = LCA(d,minimal,cutoff)
            write(g,"$(qseqid)\t$(assignment)")
            qseqid = next_qseqid
            bitscores = Dict{Taxon,Int}()
        end

        rows = next_rows
    end
    close(f)
    close(g)
end