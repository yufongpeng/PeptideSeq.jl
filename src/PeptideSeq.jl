module PeptideSeq
using IterTools,  PrettyTables

export Protein, Peptide, Fragments, digest!, modify!, ionize!, fragmentation

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
