struct BlastResult 
    qseqid::String
    sseqid::String
    pident::Float64
    length::Int
    mismatch::Int
    gapopen::Int
    qstart::Int
    qend::Int
    sstart::Int
    send::Int
    evalue::Float64
    bitscore::Float64
end


function BlastResult(line::AbstractString)
    cols = split(line, "\t")

    qseqid = cols[1]
    sseqid = cols[2]
    pident = parse(Float64, cols[3])
    length = parse(Int, cols[4])
    mismatch = parse(Int, cols[5])
    gapopen = parse(Int, cols[6])
    qstart = parse(Int, cols[7])
    qend = parse(Int, cols[8])
    sstart = parse(Int, cols[9])
    send = parse(Int, cols[10])
    evalue = parse(Float64, cols[11])
    bitscore = parse(Float64, cols[12])

    return BlastResult(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)
end

qseqid(record::BlastResult) = record.qseqid
sseqid(record::BlastResult) = record.sseqid

function blastLCA(;filepath::AbstractString, outpath::AbstractString, sqlite::SQLite.DB, taxonomy::Taxonomy.DB, method::Function, header::Bool=false, ranks=[:superkingdom, :phylum, :class, :order, :family, :genus, :species], rmselfhit=true) 
    f = open(filepath,"r")
    o = open(outpath,"w")
    try
        lca_ch = blastLCA(f; sqlite=sqlite, taxonomy=taxonomy, method=method, header=header, ranks=ranks, rmselfhit=rmselfhit)
        for (qseqid, taxon, lineage) in lca_ch
            id = taxid(taxon)
            lineage_txt = sprint(io -> print_lineage(io, lineage))
            write(o, "$(qseqid)\t$(id)\t$(lineage_txt)\n")
        end
    finally
        close(f)
        close(o)
    end
end

function blastLCA(f::IO; sqlite::SQLite.DB, taxonomy::Taxonomy.DB, method::Function, header::Bool=false, ranks=[:superkingdom, :phylum, :class, :order, :family, :genus, :species], rmselfhit=true)

    blastresult_ch = Channel{BlastResult}(500)
    lcainput_ch = Channel{Tuple{String,Dict{Taxon,BlastResult}}}(500)
    lineage_ch = Channel{Tuple{String,Taxon,Lineage}}(500)

    t1 = @async parse_blastresult!(blastresult_ch, f, header, rmselfhit)
    t2 = @async put_blastresults!(lcainput_ch, blastresult_ch, sqlite, taxonomy)
    t3 = @async lca_blastresults!(lineage_ch, lcainput_ch, method, ranks)

    bind(blastresult_ch, t1)
    bind(lcainput_ch, t2)
    bind(lineage_ch, t3)

    return lineage_ch
end

function parse_blastresult!(out_channel::Channel{BlastResult}, f::IO, header::Bool, rmselfhit::Bool)
    header ? readline(f) : nothing
    for line in eachline(f)
        record = BlastResult(line)
        if rmselfhit && qseqid(record) == sseqid(record)
            continue
        end 
        put!(out_channel, record)
    end
    close(out_channel)
end

function put_blastresults!(out_channel::Channel{Tuple{String,Dict{Taxon,BlastResult}}}, in_channel::Channel{BlastResult}, sqlite::SQLite.DB, taxonomy::Taxonomy.DB)

    results = Dict{Taxon,BlastResult}()

    while true
        record = take!(in_channel)
        taxid = get(sqlite, sseqid(record), nothing)

        if taxid === nothing
            @warn "record $(sseqid(record)) has no taxid in $(sqlite.file)"
            taxon = nothing
        else
            taxon = get(taxid, taxonomy, nothing)
        end

        if taxon === nothing
            @warn "There is no taxon correspondinig to $(taxid)!\ncontinue..."
        else
            if ! haskey(results, taxon)
                results[taxon] = record
            else
                if record.bitscore > results[taxon].bitscore
                    results[taxon] = record
                end
            end
        end

        try
            next = fetch(in_channel)
            if qseqid(record) != qseqid(next)
                put!(out_channel, (qseqid(record),results))
                results = Dict{Taxon,BlastResult}()
            end
        catch e
            put!(out_channel, (qseqid(record),results))
            break
        end
    end
    close(out_channel)
end

function lca_blastresults!(out_channel::Channel{Tuple{String,Taxon,Lineage}}, in_channel::Channel{Tuple{String,Dict{Taxon,BlastResult}}}, method::Function, ranks)
    for results in in_channel
        records = results[2]
        taxon = method(records)
        lineage = reformat(Lineage(taxon), ranks)
        put!(out_channel, (results[1], taxon, lineage))
    end
    close(out_channel)
end