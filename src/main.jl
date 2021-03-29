struct BlastResult 
    qseqid::String
    sseqid::String
    pident::Float64
    length::Int
    mismatch::Int
    gapopen::Int
    qstart::Int
    qend::Int
    sstrat::Int
    send::Int
    evalue::Float64
    bitscore::Int
end


function BlastResult(line::AbstractString)
    cols = split(line, "\t")
    return BlastResult(cols)
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
        @assert taxid !== nothing

        taxon = Taxon(taxid,taxonomy)
        if ! haskey(results, taxon)
            results[taxon] = record
        else
            if record.bitscore > results[taxon].bitscores
                results[taxon] = record
            end
        end

        next = readline(f)
        next_record = split(line, "\t")[1]

        if record.qseqid != next_qseqid  
            lineage = method(results)
            write(o,"$(qseqid)\t$(assignment)")
            results = Dict{Taxon,BlastResult}() #initialize
        end
        
        isempty(next) ? break : line = next
    end
end