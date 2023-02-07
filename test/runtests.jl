using BlastLCA
using Test
using Taxonomy
using AbstractTrees

taxonomy = Taxonomy.DB("./data/nodes.dmp", "./data/names.dmp")

@testset "parse.jl" begin
    result = BlastResult("test1\tWP_034815402.1\t100\t158\t0\t0\t1\t158\t1\t158\t3.11e-109\t318\t1207058;2630699", 1, 13, 3, 12)
    
    @test BlastLCA.qseqid(result) == "test1"
    @test BlastLCA.staxids(result) == [1207058, 2630699]
    @test BlastLCA.pident(result) == 100
    @test BlastLCA.bitscore(result) == 318

    lines = """
            test1\tWP_034815402.1\t100\t158\t0\t0\t1\t158\t1\t158\t3.11e-109\t318\t1207058;2630699
            test1\tMBL4878486.1\t96.8\t158\t5\t0\t1\t158\t1\t158\t3.49e-106\t310\t
            test1\tMAN67163.1\t96.2\t158\t6\t0\t1\t158\t1\t158\t1.00e-105\t309\t
            test1\tWP_035552119.1\t94.9\t158\t8\t0\t1\t158\t1\t158\t5.80e-105\t307\t1280948;1986606
            test1\tMBR9808867.1\t88.0\t158\t19\t0\t1\t158\t1\t158\t4.22e-98\t290\t
            test1\tWP_034825099.1\t88.0\t158\t19\t0\t1\t158\t1\t158\t8.51e-98\t289\t1280941
            test1\tWP_034794047.1\t87.3\t158\t20\t0\t1\t158\t1\t158\t1.46e-96\t286\t1280946
            """
    
    f = IOBuffer()
    write(f,lines)
    seek(f, 0)
    blastresult_ch = Channel{BlastResult}(500)
    BlastLCA.parse_blastresult!(blastresult_ch, f, false, 1, 13, 3, 12)
    
    result1 = take!(blastresult_ch)
    @test BlastLCA.qseqid(result1) == "test1"
    @test BlastLCA.staxids(result1) == [1207058, 2630699]
    @test BlastLCA.pident(result1) == 100
    @test BlastLCA.bitscore(result1) == 318

    result2 = take!(blastresult_ch)
    @test BlastLCA.qseqid(result2) == "test1"
    @test isempty(BlastLCA.staxids(result2))
    @test BlastLCA.pident(result2) == 96.8
    @test BlastLCA.bitscore(result2) == 310

    f = IOBuffer()
    write(f,lines)
    seek(f, 0)
    blastresult_ch = Channel{BlastResult}(500)
    BlastLCA.parse_blastresult!(blastresult_ch, f, false, 1, 13, 3, 12)
     
    lcainput_ch = Channel{Tuple{String,Dict{Taxon,BlastResult}}}(500)
    BlastLCA.put_blastresults!(lcainput_ch, blastresult_ch, taxonomy)
    taxon_and_result = take!(lcainput_ch)
    qid = first(taxon_and_result)
    result = last(taxon_and_result)
    @test qid == "test1"
    @test length(result) == 4
end

@testset "tree.jl" begin
    tree = BlastLCA.topolgoy(Taxon.([9593, 9605, 9606, 9597, 9601]))

    @test taxid(tree.node) == 9604
end