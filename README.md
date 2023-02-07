# BlastLCA

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://banhbio.github.io/BlastLCA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://banhbio.github.io/BlastLCA.jl/dev)
[![Build Status](https://www.travis-ci.com/banhbio/BlastLCA.jl.svg?token=TnLbrgdWxoQMPrAZynWc&branch=main)](https://travis-ci.com/banhbio/BlastLCA.jl)

This is a julia package to identify lowest common ancestor (LCA) from blast results.

It requires
- Blast result (include qseqid, bitscore, pidents and staxids)
- NCBI Taxonomy database (see [Taxonomy.jl](https://github.com/banhbio/Taxonomy.jl))

This tool is referred to [Carradec et al., 2018](https://doi.org/10.1038/s41467-017-02342-1). See Methods section (Taxonomic assignment).

## Example
```julia
using BlastLCA
using Taxonomy

db = Taxonomy.DB("/db/taxonomy/nodes.dmp", "/db/taxonomy/names.dmp")

blastout_path = "your_blast_result.txt"

minimal = 0.90
cutoff = 0.67
precision = Dict{Symbol, Float64}(
    :class => 0.50,
    :order => 0.65,
    :family => 0.80,
    :genus => 0.95,
    :species => 1.0)

ranks = [:superkingdom, :phylum, :class, :order, :family, :genus, :species]

f = x-> weightedLCA(x, minimal, cutoff, ranks, precision)

blastLCA(diamond_path, "lca_output_path.txt";
               taxonomy=db,
               method=f,
               header=false,
               qseqid_pos=1,
               staxids_pos=13,
               pident_pos=3,
               bitscore_pos=12
               );
```


## CLI

This package also support CLI.

### install

```bash
git clone https://github.com/banhbio/BlastLCA.jl.git
cd BlastLCA.jl
julia --project deps/builds.jl

export PATH="~/.julia/bin"
#or
ln -s ~/.julia/bin/blastlca /where/in/PATH
```

### usage
```bash
blastlca -h
blastlca weighted -h
```

Run weighted LCA
```bash
balstlca wighted /your/blastresult.tsv -o output.tsv --nodes-path /db/nodes.dmp --names-path /db/names.dmp
```
