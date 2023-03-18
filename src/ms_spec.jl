# Function related mass spectrometry

"""
    ionize!(protein, adducts...)

Ionize the peptides as the `adducts`.
Available adducts are "[M]", "[M+H]+", "[M+2H]2+", "[M-H]-", "[M-2H]2-".
"""
function ionize!(protein::Protein, adducts::String...)
    adduct_fn = (CONFIG["ACCURATE"] ? first : last)(ADDUCT_FN)
    peps = Vector{Peptide}(undef, length(protein.peptides) * length(adducts))
    i = 1
    for pep in protein.peptides
        for add in adducts
            peps[i] = Peptide(pep.origin, pep.position, pep.mass, add => adduct_fn[add], pep.miss_cleavage, pep.modification)
            i += 1
        end
    end
    protein.peptides = peps
    protein
end

# a: -COOH          => -CO
# b: -OH            => 
# c: -H2O + NH4     => +NH3
# x: -H2O + COOH    => +CO - 2H + H2O   / peptide - a - 2H
# y: +H             => +H2O             / peptide - b
# z: -NH2           => -NH3 + H2O       / peptide - c

"""
    fragmentation(peptide; ion_type = [:b, :y], charge_state = :auto)

Fragmentation of a peptide.
`ion_type` is the type of fragments. It can be an vector containing :a, :b, :c, :x, :y, :z. The default is [:b, :c] which are major fragments in CID or HCD.
`charge_state` determines the number of charges on the fragments. It can be a vector containing integers or a symbol. The default is `:auto` which means that doulbly charged fragments will be included for fragments containing more than 5 amino acids.
"""
function fragmentation(peptide::Peptide; ion_type = [:b, :y], charge_state = :auto)
    fn = CONFIG["ACCURATE"] ? first : last
    aa_ms = fn(AA_MS)
    add_ms = fn(ADD_MS)
    modification_ms = fn(MODIFICATION_MS)
    adduct_fn = fn(ADDUCT_FN)

    neutral_fragments = Dict{Symbol, Vector{Float64}}()
    seq_mass = _modified_mass(peptide, aa_ms, modification_ms)
    if :a in ion_type
        a = a_ion(seq_mass, add_ms)
        push!(neutral_fragments, :a => a)
        :x in ion_type && push!(neutral_fragments, :x => reverse(peptide.mass .- a .- h_ms * 2))
    elseif :x in ion_type
        push!(neutral_fragments, :x => x_ion(seq_mass, add_ms))
    end
    if :b in ion_type
        b = b_ion(seq_mass, add_ms)
        push!(neutral_fragments, :b => b)
        :y in ion_type && push!(neutral_fragments, :y => reverse(peptide.mass .- b))
    elseif :y in ion_type
        push!(neutral_fragments, :y => y_ion(seq_mass, add_ms))
    end
    if :c in ion_type
        c = c_ion(seq_mass, add_ms)
        push!(neutral_fragments, :c => c)
        :z in ion_type && push!(neutral_fragments, :z => reverse(peptide.mass .- c))
    elseif :z in ion_type
        push!(neutral_fragments, :z => z_ion(seq_mass, add_ms))
    end

    _, precursor_charge_state = match(r"(\[.*\])(.*)", first(peptide.adduct))
    Fragments(peptide, last(peptide.adduct)(peptide.mass), _charge_fragments(neutral_fragments, charge_state, precursor_charge_state, adduct_fn))
end

function _charge_fragments(neutral_fragments::Dict{Symbol, Vector{Float64}}, charge_state, precursor_charge_state, adduct_fn)
    charge, ion_mode = match(r"(\d*)([+, -])", precursor_charge_state)
    maxcharge = charge == "" ? 1 : parse(Int, charge)
    maximum(collect(charge_state)) > maxcharge && throw(ArgumentError("Fragments can not have more charges than precursor; try differnt `charge_state"))
    ion_type = collect(keys(neutral_fragments))
    sort!(ion_type, by = x -> findfirst(==(x), [:a, :x, :b, :y, :c, :z]))
    n = length(ion_type) * length(neutral_fragments[ion_type[1]]) * length(charge_state)
    mass = Vector{Float64}(undef, n)
    type = Vector{String}(undef, n)
    id = 1
    for ion in ion_type 
        for ncharge in charge_state 
            charge = ncharge == 1 ? "" : "$ncharge"
            adduct = "[M$(ion_mode)$(charge)H]$(charge)$(ion_mode)"
            charge = UPPER_INDEX[charge * ion_mode]
            for (i, v) in enumerate(neutral_fragments[ion])
                type[id] = String(ion) * Char(8320 + i) * charge
                mass[id] = adduct_fn[adduct](v)
                id += 1
            end
        end
    end
    (type = type, mass = mass)
end

function _charge_fragments(neutral_fragments::Dict{Symbol, Vector{Float64}}, charge_state::Symbol, precursor_charge_state, adduct_fn)
    scharge, ion_mode = match(r"(\d*)([+, -])", precursor_charge_state)
    maxcharge = scharge == "" ? 1 : parse(Int, scharge)
    ion_type = collect(keys(neutral_fragments))
    sort!(ion_type, by = x -> findfirst(==(x), [:a, :x, :b, :y, :c, :z]))
    n_seq = length(neutral_fragments[ion_type[1]]) 
    n = length(ion_type) * (n_seq + max(n_seq - 5, 0) * (maxcharge - 1))
    mass = Vector{Float64}(undef, n)
    type = Vector{String}(undef, n)
    charge_state == :auto || return
    id = 1
    for ion in ion_type
        charges = [""]
        for (i, v) in enumerate(neutral_fragments[ion])
            for charge in charges
                adduct = "[M$(ion_mode)$(charge)H]$(charge)$(ion_mode)"
                charge = UPPER_INDEX[charge * ion_mode]
                type[id] = String(ion) * Char(8320 + i) * charge
                mass[id] = adduct_fn[adduct](v)
                id += 1
            end
            i == 5 && maxcharge > 1 && push!(charges, "2")
        end
    end
    (type = type, mass = mass)
end

function _modified_mass(peptide::Peptide, aa_ms, modification_ms)
    seq_mass = map(collect(peptide.origin[peptide.position])) do c
        aa_ms[c]
    end
    for (k, v) in peptide.modification
        id = filter(in(v), peptide.position)
        isempty(id) && continue
        seq_mass[id] .+= modification_ms[k]
    end
    seq_mass
end

a_ion(seq_mass, add_ms) = accumulate(+, seq_mass)[1:end - 1] .- add_ms["CO"]
b_ion(seq_mass, add_ms) = accumulate(+, seq_mass)[1:end - 1]
c_ion(seq_mass, add_ms) = accumulate(+, seq_mass)[1:end - 1] .+ add_ms["NH3"]
x_ion(seq_mass, add_ms) = @. accumulate(+, reverse(seq_mass))[1:end - 1] + add_ms["CO"] - 2 * h_ms + add_ms["DI"]
y_ion(seq_mass, add_ms) = accumulate(+, reverse(seq_mass))[1:end - 1] .+ add_ms["DI"]
z_ion(seq_mass, add_ms) = @. accumulate(+, seq_mass)[1:end - 1] - add_ms["NH3"] + add_ms["DI"]