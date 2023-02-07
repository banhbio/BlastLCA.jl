using Comonicon
using JSON

"""
calcurate LCA using weighted LCA method.

# Args

- `input`: balast output file. should be in the tsv format.

# Options

- `-o, --output <path>`: is the path to the output file. `Required`
- `--nodes-path <path>`: is the path to the nodes file. `Required`
- `--names-path <path>`: is the path to the names file. `Required`
- `--minimal <float>`: 
- `--cutoff <float>` :
- `--ranks <ranks>`: should be in the JSON format. Default `"[\"superkingdom\", \"phylum\", \"class\", \"order\", \"family\", \"genus\", \"species\"]"`
- `--precision <precision>`: should be in the JSON format. Default `"{\"class\" : 0.50, \"order\" : 0.65, \"family\" : 0.80, \"genus\" : 0.95, \"species\" : 1.0}"`
- `--qseqid-pos <int>`: Default `1`
- `--staxids-pos <int>`: Default `13`
- `--pident-pos <int>`: Default `3`
- `--bitscore-pos <int>`: Default `12`
"""
@cast function weighted(input;
                output=nothing,
                nodes_path=nothing,
                names_path=nothing,
                minimal::Float64 = 0.90,
                cutoff::Float64 = 0.67,
                ranks="[\"superkingdom\", \"phylum\", \"class\", \"order\", \"family\", \"genus\", \"species\"]",
                precision ="{\"class\" : 0.50, \"order\" : 0.65, \"family\" : 0.80, \"genus\" : 0.95, \"species\" : 1.0}",
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

        f = x -> weightedLCA(x, minimal, cutoff, r, pre)

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
        cmd_error("`$(e.msg)`")
    end
end

"""
calcurate lowest common ancestors from blast results. see https://github.com/banhbio/BlastLCA.jl
"""
@main
