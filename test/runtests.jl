using PeptideSeq
using Test

@testset "PeptideSeq.jl" begin
    s = "DPCHKPKRRKP"
    p1 = Protein(s)
    p2 = Protein(s)
    digest!(p1, 2, "Trypsin")
    modify!(p1, "3NPH")
    modify!(p2, "3NPH")
    digest!(p2, 2, "Trypsin")
    ionize!(p1, "[M+H]+", "[M+2H]2+")
    ionize!(p2, "[M+H]+", "[M+2H]2+")
    fragmentation(p1.peptides[1])
    fragmentation(p2.peptides[1])
    fragmentation(p1.peptides[2])
    fragmentation(p2.peptides[2])
end
