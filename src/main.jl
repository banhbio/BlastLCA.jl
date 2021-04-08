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

function blastLCA(;filepath::AbstractString, outpath::AbstractString, sqlite::SQLite.DB, taxonomy::Taxonomy.DB, method::Function, header::Bool=false)
    f = open(filepath,"r")
    o = open(outpath,"w")
    blastLCA(f, o; sqlite=sqlite, taxonomy=taxonomy, method=method, header=header)
    close(f)
    close(o)
end

function blastLCA(f::IOStream, o::IOStream; sqlite::SQLite.DB, taxonomy::Taxonomy.DB, method::Function, header::Bool=false)
    header ? readline(f) : nothing #skip header

    results = Dict{Taxon,BlastResult}()

    line = readline(f) #read first line
    while true
        record = BlastResult(line)
        taxid = get(sqlite, record.sseqid, nothing)
        if taxid === nothing
            @warn "record $(record.sseqid) has no taxid in $(sqlite.file)"
            taxon = nothing
        else
            taxon = get(taxid, taxonomy, nothing)
        end


        if taxon === nothing
            @warn "There is no taxon correspondinig to $(taxid)!"
            @warn "continue..."
        else
            if ! haskey(results, taxon)
                results[taxon] = record
            else
                if record.bitscore > results[taxon].bitscore
                    results[taxon] = record
                end
            end
        end

        next = readline(f)
        next_qseqid = split(next, "\t")[1]
        if record.qseqid != next_qseqid
            lineage = method(results)
            lineage_txt = lineage_line(lineage)
            write(o,"$(record.qseqid)\t$(lineage_txt)\n")
            results = Dict{Taxon,BlastResult}() #initialize
        end
        
        isempty(next) ? break : line = next
    end
end