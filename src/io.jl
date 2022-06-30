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