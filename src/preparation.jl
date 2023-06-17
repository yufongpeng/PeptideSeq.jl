# Digestion and modification

"""
    modify!(protein, modification...)

Modify an intact protein or digested protein.
Modification must be a key of `MODIFICATION_SITE` or it will be ignored.
See object `MODIFICATION_SITE` for available modifications.
Currently, "3NPH" is supported.
"""
function modify!(protein::Protein, modification::String...)
    for k in modification
        haskey(MODIFICATION_SITE, k) || continue
        locs = Int[]
        for loc in MODIFICATION_SITE[k]
            if isa(loc, Regex)
                append!(locs, locc.offset for locc in eachmatch(loc, protein.origin))
            elseif loc == "^"
                push!(locs, 0)
                push!(locs, 1)
            elseif loc == "\$"
                push!(locs, length(protein.origin))
                push!(locs, -1)
            else
                append!(locs, findall(==(first(loc)), protein.origin))
            end
        end
        protein.modification[k] = locs
    end
    # If digestion had been done, add mass to each peptides
    return isempty(protein.peptides) ? protein : _modify_mass!(protein, modification...)
end

function modify!(protein::Protein, modification::Dict{String, Vector{Int}})
    for (k, v) in modification
        protein.modification[k] = sort!(union!(get(protein.modification, k, Int[]), v))
    end
    protein
end

function _modify_mass!(protein::Protein, modification::String...)
    modification_ms = (CONFIG["ACCURATE"] ? first : last)(MODIFICATION_MS)

    for (k, v) in protein.modification
        for (i, pep) in enumerate(protein.peptides)
            id = filter(in(v), pep.position)
            in(0, v) && !=(1, first(pep.position)) && pushfirst!(id, first(pep.position))
            in(-1, v) && !=(length(protein.peptides), last(pep.position)) && push!(id, last(pep.position))
            isempty(id) && continue
            protein.peptides[i].modification[k] = sort!(union!(get(protein.peptides[i].modification, k, Int[]), id))
            n = length(protein.peptides[i].modification[k])
            protein.peptides[i].mass += n * modification_ms[k]
        end
    end
    protein
end

function full_digestion(sequence::AbstractString, enzyme::String)
    regex = ENZYME[enzyme]
    cleavage_sites = mapreduce(union, regex) do rs
        [s.offset for s in eachmatch(rs, sequence, overlap = true)]
    end
    pushfirst!(cleavage_sites, 0)
    push!(cleavage_sites, length(sequence))
    map(@views zip(cleavage_sites[1:end - 1], cleavage_sites[2:end])) do (s, e)
        (s + 1):e
    end
end

function merge_modification(modifications)
    new = Dict{String, Vector{Int}}()
    for d in modifications
        for (k, v) in d
            new[k] = union(get(new, k, Int[]), v)
        end
    end
    for v in values(new)
        sort!(v)
    end
    new
end

"""
    digest!(protein, n_miss, enzyme = "")

Digest a protein with an enzyme.
`n_miss` is the allowed number of miss cleavage.
If `protein.enzyme` is an empty string, `enzyme` must be provided.
See object `ENZYME` for available enzymes.
"""
function digest!(protein::Protein, n_miss::Int, enzyme::String = "")
    fn = CONFIG["ACCURATE"] ? first : last
    aa_ms = fn(AA_MS)
    di_ms = fn(ADD_MS)["DI"]
    modification_ms = fn(MODIFICATION_MS)
    isempty(enzyme) ? (isempty(protein.enzyme) && throw(ArgumentError("Please provide enzyme!"))) : (protein.enzyme = enzyme)

    seq_id = full_digestion(protein.origin, protein.enzyme)
    full_digestion_mass = map(seq_id) do id
        mapreduce(+, protein.origin[id]) do c
            aa_ms[c]
        end
    end    
    adduct = "[M]" => identity
    n_seq = length(seq_id)
    n_miss += 1
    n_pep = round(Int, (2 * n_seq - n_miss + 1) * n_miss / 2) 
    mods = Dict{String, Vector{Int}}()
    modification = map(1:n_seq) do _
        deepcopy(mods)
    end

    # If modification had been done, add mass to each peptides
    # split out ?
    if !isempty(protein.modification)
        for (k, v) in protein.modification
            for (i, seq) in enumerate(seq_id)
                id = filter(in(v), seq)
                isempty(id) && continue
                push!(modification[i], k => id)
                n = length(id)
                full_digestion_mass[i] += n * modification_ms[k]
            end
        end
    end
    
    protein.peptides = Vector{Peptide}(undef, n_pep)
    protein.peptides[1:n_seq] .= Peptide.(protein.origin, seq_id, full_digestion_mass .+ di_ms, adduct, 0, modification)

    n_miss > 1 || return protein
    i = n_seq + 1
    for imer in 2:n_miss
        prev_pep = reduce(+, full_digestion_mass[1:imer - 1]) + di_ms
        prev_aa = 0
        for (id, ids) in zip(imer:n_seq, IterTools.partition(1:n_seq, imer, 1))
            prev_pep += full_digestion_mass[id] - prev_aa
            protein.peptides[i] = Peptide(  protein.origin, 
                                            seq_id[first(ids)].start:seq_id[last(ids)].stop, 
                                            prev_pep, 
                                            adduct, 
                                            imer - 1, 
                                            merge_modification(modification[first(ids):last(ids)])
                                            )
            prev_aa = full_digestion_mass[id - imer + 1]
            i += 1
        end
    end
    protein
end