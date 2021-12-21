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

function blastLCA(;filepath::AbstractString, outpath::AbstractString, sqlite::SQLite.DB, taxonomy::Taxonomy.DB, method::Function, header::Bool=false, ranks=[:superkingdom, :phylum, :class, :order, :family, :genus, :species]) 
    f = open(filepath,"r")
    o = open(outpath,"w")
    try
        blastLCA(f, o; sqlite=sqlite, taxonomy=taxonomy, method=method, header=header, ranks=ranks)
    finally
        close(f)
        close(o)
    end
end

function blastLCA(f::IOStream, o::IOStream; sqlite::SQLite.DB, taxonomy::Taxonomy.DB, method::Function, header::Bool=false, ranks=[:superkingdom, :phylum, :class, :order, :family, :genus, :species])

    blastresult_ch = Channel{BlastResult}(32)
    lcainput_ch = Channel{Tuple{String,Dict{Taxon,BlastResult}}}(32)
    lineage_ch = Channel{Tuple{String,Lineage}}(32)

    @async parse_blastresult!(blastresult_ch, f, header)
    @async put_blastresults!(lcainput_ch, blastresult_ch, sqlite, taxonomy)
    @async lca_blastresults!(lineage_ch, lcainput_ch, method, ranks)

    for (qseqid, lineage) in lineage_ch
        lineage_txt = sprint(io -> print_lineage(io, lineage))
        write(o, "$(qseqid)\t$(lineage_txt)\n")
    end
end

function parse_blastresult!(out_channel::Channel{BlastResult}, f::IOStream, header::Bool)
    header ? readline(f) : nothing
    for line in eachline(f)
        record = BlastResult(line)
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
                results = Dict{Taxon, BlastResult}()
            end
        catch e
            put!(out_channel, (qseqid(record),results))
            break
        end
    end
    close(out_channel)
end

function lca_blastresults!(out_channel::Channel{Tuple{String, Lineage}}, in_channel::Channel{Tuple{String, Dict{Taxon, BlastResult}}}, method::Function, ranks)
    for results in in_channel
        taxon = method(results[2])
        lineage = reformat(Lineage(taxon), ranks)
        put!(out_channel, (results[1], lineage))
    end
    close(out_channel)
end