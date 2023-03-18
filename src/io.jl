# Custom enzyme and modification
"""
    add_enzyme!(source)

Add custom enzyme. The `source` must be a tsv file. The first row is the header (can be empty), the first column is the name of enzyme and the other column is regular expressions(r"...") of the cleavage sites. 
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

Add custom modification. The `source` must be a tsv file. The first row is the header (can be empty), the first column is the name of modification, the second and third columns are addtional monoisotopic mass and average mass, respectively, and the other columns are the modification site. 
Modification sites can be a string without quotation or regular expression like r"...". ^ repressents the N-terminal and \$ repressents the C-terminal. 
See the example file "config_example/MODIFICATION.tsv".
"""
function add_modification!(source)
    for m in CSV.Rows(source, delim = "\t")
        push!(MODIFICATION_MS[1], m[1] => parse(Float64, m[2]))
        push!(MODIFICATION_MS[2], m[1] => parse(Float64, m[3]))
        push!(MODIFICATION_SITE, m[1] => [startwith(loc, "r") ? eval(Meta.parse(loc)) : loc for loc in getindex.(Ref(m), 4:length(m)) if !ismissing(loc)])
    end
end

function read_msdial(source)
    return
end

# Custom display
function display_modification(modification; position::Bool = false)
    s = position ? [k * " (@ " * join(v, ", ") * ")" for (k, v) in modification] : keys(modification)
    join(s, ", ")
end

function print_adduct(adduct::String)
    body, charge = match(r"(\[.*\])(.*)", adduct)
    if isempty(charge)
        return body
    end
    charge = UPPER_INDEX[charge]
    return body * charge
end

import Base: show

function Base.show(io::IO, protein::Protein)
    isempty(protein.peptides) && return 
    dt = [( Peptide_sequence    = protein.origin[pep.position], 
            Position            = pep.position,
            Mass                = last(pep.adduct)(pep.mass),
            Adduct              = print_adduct(first(pep.adduct)),
            Miss_cleavage       = pep.miss_cleavage,
            Modification        = display_modification(pep.modification; position = true)) for pep in protein.peptides]
    pretty_table(io, dt, title = "Peptides: ", 
                header = collect(propertynames(dt[1])), 
                header_alignment = :L, 
                alignment=[:r, :r, :r, :l, :r, :l])
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, protein::Protein) 
    printstyled(io, rpad("Protein: ", 14); bold = true)
    print(io, protein.origin, "\n")
    printstyled(io, "Modification: "; bold = true)
    print(io, display_modification(protein.modification), "\n")
    printstyled(io, rpad("Enzyme: ", 14); bold = true)
    print(io, protein.enzyme, "\n")
    show(io, protein)
end

function Base.show(io::IO, fragments::Fragments)
    pretty_table(io, fragments.fragments, title = "Fragment ions: ", 
                    header = [:Type, :Mass], 
                    header_alignment = :L, 
                    alignment=[:l, :r])
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, fragments::Fragments) 
    pep = fragments.peptide
    printstyled(io, rpad("Peptide: ", 14); bold = true)
    print(io, pep.origin[pep.position], "\n")
    printstyled(io, "Modification: "; bold = true)
    print(io, display_modification(pep.modification; position = true), "\n")
    printstyled(io, rpad("Precursor: ", 14); bold = true)
    print(io, fragments.precursor, "\n")
    printstyled(io, rpad("Adduct: ", 14); bold = true)
    print(io, print_adduct(first(pep.adduct)), "\n")
    show(io, fragments)
end