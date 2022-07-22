module PeptideSeq
using IterTools,  PrettyTables

export Protein, Peptide, Fragments, 

    digest!, modify!, ionize!, fragmentation,

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

"""
    add_enzyme!(source)

Add custom enzyme. The `source` must be a tsv file. The first row is the header (can be empty), the first column is the name of enzyme and the other column is regular expressions of the cleavage sites. 
If multiple regular expressions are given, the digested sites will be the union of all possible sites.
See the example file "config_example/ENZYME.tsv".
"""
function add_enzyme!(source)
    for e in CSV.Rows(source, delim = "\t")
        push!(ENZYME,  e[1] => [eval(Meta.parse(r)) for r in getindex.(Ref(e), 2:length(e)) if !ismissing(r)])
    end
end

"""
    add_modification!(source)

Add custom modification. The first row is the header (can be empty), the first column is the name of modification, the second and third columns are addtional monoisotopic mass and average mass, respectively, and the other columns are the modification site. 
Modification sites can be string without quotation or regular expression like r"...". ^ repressents the N-terminal and \$ repressents the C-terminal. 
See the example file "config_example/MODIFICATION.tsv".
"""
function add_modification!(source)
    for m in CSV.Rows(source, delim = "\t")
        push!(MODIFICATION_MS[1], m[1] => parse(Float64, m[2]))
        push!(MODIFICATION_MS[2], m[1] => parse(Float64, m[3]))
        push!(MODIFICATION_SITE, m[1] => [startwith(loc, "r") ? eval(Meta.parse(loc)) : loc for loc in getindex.(Ref(m), 4:length(m)) if !ismissing(loc)])
    end
end


end
