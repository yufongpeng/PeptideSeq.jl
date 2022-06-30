# Function related mass spectrometry

function ionize!(protein::Protein, adducts::String...)
    if CONFIG["ACCURATE"]
        adduct_fn = first(ADDUCT_FN)
    else
        adduct_fn = last(ADDUCT_FN_AVERAGE)
    end
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

function fragmentation(peptide::Peptide; ion_type = [:b, :y], charge_state = :auto)
    if CONFIG["ACCURATE"]
        aa_ms = first(AA_MS)
        add_ms = first(ADD_MS)
        modification_ms = first(MODIFICATION_MS)
        adduct_fn = first(ADDUCT_FN)
    else
        aa_ms = last(AA_MS)
        add_ms = last(ADD_MS)
        modification_ms = last(MODIFICATION_MS)
        adduct_fn = last(ADDUCT_FN_AVERAGE)
    end

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
    charge = charge == "" ? 1 : parse(Int, charge)
    @assert maximum(collect(charge_state)) <= charge "Fragments can not have more charges than precursor"
    ion_type = collect(keys(neutral_fragments))
    sort!(ion_type, by = x -> findfirst(==(x), [:a, :x, :b, :y, :c, :z]))
    n = length(ion_type) * length(neutral_fragments[ion_type[1]]) * length(charge_state)
    mass = Vector{Float64}(undef, n)
    type = Vector{String}(undef, n)
    id = 1
    for ion in ion_type 
        for charge in charge_state 
            for (i, v) in enumerate(neutral_fragments[ion])
                charge = charge == 1 ? "" : "$charge"
                adduct = "[M$(ion_mode)$(charge)H]$(charge)$(ion_mode)"
                charge = UPPER_INDEX[charge * ion_mode]
                type[id] = String(ion) * Char(8320 + i) * charge
                mass[id] = adduct_fn[adduct](v)
                id += 1
            end
        end
    end
    (type = type, mass = mass)
end

function _charge_fragments(neutral_fragments::Dict{Symbol, Vector{Float64}}, charge_state::Symbol, precursor_charge_state, adduct_fn)
    charge, ion_mode = match(r"(\d*)([+, -])", precursor_charge_state)
    charge = charge == "" ? 1 : parse(Int, charge)
    ion_type = collect(keys(neutral_fragments))
    sort!(ion_type, by = x -> findfirst(==(x), [:a, :x, :b, :y, :c, :z]))
    n_seq = length(neutral_fragments[ion_type[1]]) 
    n = length(ion_type) * (n_seq + max(n_seq - 5, 0) * (charge - 1))
    mass = Vector{Float64}(undef, n)
    type = Vector{String}(undef, n)
    if charge_state != :auto
        return
    end
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
            i == 5 && (charge - 1) > 0 && push!(charges, "2")
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

function a_ion(seq_mass, add_ms)
    accumulate(+, seq_mass)[1:end - 1] .- add_ms["CO"]
end

function b_ion(seq_mass, add_ms)
    accumulate(+, seq_mass)[1:end - 1]
end

function c_ion(seq_mass, add_ms)
    accumulate(+, seq_mass)[1:end - 1] .+ add_ms["NH3"]
end

function x_ion(seq_mass, add_ms)
    @. accumulate(+, reverse(seq_mass))[1:end - 1] + add_ms["CO"] - 2 * h_ms + add_ms["DI"]
end

function y_ion(seq_mass, add_ms)
    accumulate(+, reverse(seq_mass))[1:end - 1] .+ add_ms["DI"]
end

function z_ion(seq_mass, add_ms)
    @. accumulate(+, seq_mass)[1:end - 1] - add_ms["NH3"] + add_ms["DI"]
end