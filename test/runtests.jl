using PeptideSeq
using Test

@testset "PeptideSeq.jl" begin
    s = "DPCHKPKRRKP"
    p1 = Protein(s)
    p2 = Protein(s)
    digest!(p1, 2, "Trypsin")
    @test length(p1.peptides) == 9
    @test p1.peptides[1].position == 1:7
    modify!(p1, "3NPH")
    @test p1.modification["3NPH"] == [1, 11, -1]
    @test p1.peptides[1].modification["3NPH"] == [1, 7]
    modify!(p2, "3NPH")
    digest!(p2, 2, "Trypsin")
    @test p2.modification["3NPH"] == [1, 11, -1]
    @test p2.peptides[1].modification["3NPH"] == [1]
    ionize!(p1, "[M+H]+", "[M+2H]2+")
    ionize!(p2, "[M+H]+", "[M+2H]2+")
    fragmentation(p1.peptides[1])
    fragmentation(p2.peptides[1])
    fragmentation(p1.peptides[2]; charge_state = [1])
    fragmentation(p2.peptides[2]; charge_state = [1, 2])
    db = [fragmentation(p2.peptides[16]), fragmentation(p2.peptides[10])]
    data = read_msdial("data.txt")
    result = find_peptides(data, db)
    @test findall(result.var"Peptide Matched") == [121, 155]
    @test result.Peptide[121] == db[2]
    @test result.Score[121] > 0.8
    @test result.Peptide[155] == db[1]
    @test result.Score[155] < 0.1
end
