using Comonicon
using JSON

"""
"""
@main function blastlca(input;
                output=nothing,
                nodes_path=nothing,
                names_path=nothing,
                minimal::Float64 = 0.90,
                cutoff::Float64 = 0.67,
                ranks="[\"superkingdom\", \"phylum\", \"class\", \"order\", \"family\", \"genus\", \"species\"]",
                precision ="{\"class\" : 0.50, \"order\" : 0.65, \"family\" : 0.80, \"genus\" : 0.95, \"species\" : 1.0}",
                fun="weightedLCA",
                header::Bool=false,
                qseqid_pos::Int=1,
                staxids_pos::Int=13,
                pident_pos::Int=3,
                bitscore_pos::Int=12
                )

    try    
        db = Taxonomy.DB(nodes_path, names_path)

        r = [Symbol(class) for class in JSON.parse(ranks)]
        pre = Dict(Symbol(k) => v for (k, v) in Dict(JSON.parse(precision)))

        if fun=="weightedLCA"
            f = x -> weightedLCA(x, minimal, cutoff, r, pre)
        else
            error("not supported function")
        end

        blastLCA(input, output;
                    taxonomy=db,
                    method=f,
                    header=header,
                    qseqid_pos=qseqid_pos,
                    staxids_pos=staxids_pos,
                    pident_pos=pident_pos,
                    bitscore_pos=bitscore_pos
        );

    catch e
        cmd_error(e.msg)
    end
end