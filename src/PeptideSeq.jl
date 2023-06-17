module PeptideSeq
using IterTools,  PrettyTables, TypedTables

export Protein, Peptide, Fragments, 

    digest!, modify!, ionize!, fragmentation, find_peptides, read_msdial,

    add_enzyme!, add_modification!,

    MODIFICATION_SITE, ENZYME, CONFIG

"""
    Peptide

Repressentation of a peptide digested from a `Protein`.
# Field
* `origin`: Oringinal protein sequence
* `position`: Position of this peptide counting from N-terminal
* `mass`: Monoisotopic or average mass depending on `CONFIG["ACCURACY"]`
* `adduct`: Adduct of the ionized peptide. When the peptide is not ionized, it is "[M]". 
* `miss_cleavage`: Number of miss cleavages allowed
* `modification`: modification on the peptide
"""
mutable struct Peptide
    origin::AbstractString
    position::UnitRange{Int}
    mass::Float64
    adduct::Pair{String, Function}
    miss_cleavage::Int
    modification::Dict{String, Vector{Int}}
end

struct Fragments
    peptide::Peptide
    precursor::Float64
    fragments::NamedTuple{(:type, :mass), Tuple{Vector{String}, Vector{Float64}}}
end

"""
    Protein
    Protein(sequence, enzyme = "", modification = Dict{String, Vector{Int}}())

Repressentation of a protein
# Field
* `origin`: Protein sequence
* `peptides`: a vector of `Peptide` digested with `enzyme`
* `modification`: modification on the protein
* `enzyme`: an enzyme for digesting the protein
"""
mutable struct Protein
    origin::AbstractString
    peptides::Vector{Peptide}
    modification::Dict{String, Vector{Int}}
    enzyme::String
end

Protein(sequence::AbstractString, enzyme::String = "", modification::Dict{String, Vector{Int}} = Dict{String, Vector{Int}}()) = 
    Protein(sequence, Peptide[], modification, enzyme)


include("config.jl")
include("io.jl")
include("preparation.jl")
include("ms_spec.jl")

end
