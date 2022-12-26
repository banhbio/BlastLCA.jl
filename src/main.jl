struct BlastResult 
    qseqid::String
    staxids::Vector{Int}
    pident::Float64
    bitscore::Float64
end


function BlastResult(line::AbstractString, qseqid_pos::Int, staxid_pos::Int, pident_pos::Int, bitscore_pos::Int)
    cols = split(line, "\t")

    qseqid = cols[qseqid_pos]
    staxids = parse.(Int, split(cols[staxid_pos], ";"))
    pident = parse(Float64, cols[pident_pos])
    bitscore = parse(Float64, cols[bitscore_pos])

    return BlastResult(qseqid, staxids, pident, bitscore)
end

qseqid(record::BlastResult) = record.qseqid
staxids(record::BlastResult) = record.staxids
pident(record::BlastResult) = record.pident
bitscore(record::BlastResult) = record.bitscore

function blastLCA(filepath::AbstractString, outpath::AbstractString; kwargs...) 
    f = open(filepath,"r")
    o = open(outpath,"w")
    try
        lca_ch = blastLCA(f; kwargs...)
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

function blastLCA(df::DataFrame; kwargs...)
    f = IOBuffer()
    CSV.write(f, df; delim="\t")
    seek(f, 0)
    lca_ch = blastLCA(f; kwargs...)
    lca_rows = NamedTuple{(:qseqid, :taxon, :lineage), Tuple{String, Taxon, Lineage}}[]
    for (qseqid, taxon, lineage) in lca_ch
        row = (qseqid = qseqid, taxon = taxon, lineage = lineage)
        push!(lca_rows, row)
    end
    return DataFrame(lca_rows)
end

function blastLCA(f::IO; taxonomy::Taxonomy.DB, method::Function, header::Bool=false, qseqid_pos::Int=1, staxid_pos::Int=2, pident_pos::Int=3, bitscore_pos::Int=4, ranks=[:superkingdom, :phylum, :class, :order, :family, :genus, :species])

    blastresult_ch = Channel{BlastResult}(500)
    lcainput_ch = Channel{Tuple{String,Dict{Taxon,BlastResult}}}(500)
    lineage_ch = Channel{Tuple{String,Taxon,Lineage}}(500)

    t1 = @async parse_blastresult!(blastresult_ch, f, header, qseqid_pos, staxid_pos, pident_pos, bitscore_pos)
    t2 = @async put_blastresults!(lcainput_ch, blastresult_ch, taxonomy)
    t3 = @async lca_blastresults!(lineage_ch, lcainput_ch, method, ranks)

    bind(blastresult_ch, t1)
    bind(lcainput_ch, t2)
    bind(lineage_ch, t3)

    return lineage_ch
end

function parse_blastresult!(out_channel::Channel{BlastResult}, f::IO, header::Bool, qseqid_pos::Int, staxid_pos::Int, pident_pos::Int, bitscore_pos::Int)
    header ? readline(f) : nothing
    for line in eachline(f)
        record = BlastResult(line, qseqid_pos, staxid_pos, pident_pos, bitscore_pos)
        put!(out_channel, record)
    end
    close(out_channel)
end

function put_blastresults!(out_channel::Channel{Tuple{String,Dict{Taxon,BlastResult}}}, in_channel::Channel{BlastResult}, taxonomy::Taxonomy.DB)

    results = Dict{Taxon,BlastResult}()

    while true
        record = take!(in_channel)
        taxid = staxids(record) |> first
        taxon = get(taxonomy, taxid, nothing)

        if taxon === nothing
            @warn "There is no taxon correspondinig to $(taxid)!\ncontinue..."
        elseif isdescendant(taxon, Taxon(28384, taxonomy)) || isdescendant(taxon, Taxon(12908, taxonomy))
            @warn "This sequence ($(taxid)) is coming from other sequences or unclassifiedsequences\ncontinue..."
        else
            if ! haskey(results, taxon)
                results[taxon] = record
            else
                if bitscore(record) > bitscore(results[taxon])
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
        taxon = method(results)
        lineage = reformat(Lineage(taxon), ranks)
        put!(out_channel, (results[1], taxon, lineage))
    end
    close(out_channel)
end