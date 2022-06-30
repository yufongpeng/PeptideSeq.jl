# Mass and other config
using CSV

const CONFIG = Dict(
    "ACCURATE" => true
)

const AA_MS = (
    Dict(
        'A' => 71.03711,
        'R' => 156.10111,
        'N' => 114.04293,
        'D' => 115.02694,
        'C' => 103.00919,
        'E' => 129.04259,
        'Q' => 128.05858,
        'G' => 57.02146,
        'H' => 137.05891,
        'I' => 113.08406,
        'L' => 113.08406,
        'K' => 128.09496,
        'M' => 131.04049,
        'F' => 147.06841,
        'P' => 97.05276,
        'S' => 87.03203,
        'T' => 101.04768,
        'W' => 186.07931,
        'Y' => 163.06333,
        'V' => 99.06841,
    ), 
    Dict(
        'A' => 71.0788,
        'R' => 156.1875,
        'N' => 114.1038,
        'D' => 115.0886,
        'C' => 103.1388,
        'E' => 129.1155,
        'Q' => 128.1307,
        'G' => 57.0519,
        'H' => 137.1411,
        'I' => 113.1594,
        'L' => 113.1594,
        'K' => 128.1741,
        'M' => 131.1926,
        'F' => 147.1766,
        'P' => 97.1167,
        'S' => 87.0782,
        'T' => 101.1051,
        'W' => 186.2132,
        'Y' => 163.1760,
        'V' => 99.1326
    )
)
const ADD_MS = (
    Dict(
        "DI" => 18.010565,
        "H" => 1.007825,
        "CO" => 27.994915,
        "NH3" => 17.026549
    ),
    Dict(
        "DI" => 18.01528,
        "H" => (1.00784+1.00811) / 2,
        "CO" => 27.994915,
        "NH3" => 17.026549
    )  
)

const MODIFICATION_SITE = Dict{String, Vector{String}}("3NPH" => ["D", "E", "\$"])

const MODIFICATION_MS = (
    Dict{String, Float64}("3NPH" => 135.043262),
    Dict{String, Float64}("3NPH" => 135.12472)
)

const ENZYME = Dict(
    "Trypsin" => [r"[K, R][^P]"],                           # R/KX, X≠P
    "Trypsin(low specificity)" =>  [r"[K, R]."],            # R/KX
    "Trypsin(high specificity)" => [r"K[^P]", r"R[^P, R]"], # R/KX, X≠P, no RR
    "Trypsin/CNBr" => [r"[K, R, M][^P]"],                   # R/K/MX, X≠P
    "Lys C" => [r"K."], # KX
    "Lys N" => [r".K"], # XK
    "CNBr" => [r"M."],  # MX
    "Arg C" => [r"R[^P]"], # RX, X≠P
    "Asp N" => [r".D"], # XD
    "Glu N" => [r".E"], # XE
    "Asp N/Glu N" => [r".D", r".E"],
    "Asp N/Lys C" => [r".D", r"K."],
    "Asp N/Glu N/Lys C" => [r".D", r".E", r"K."],
    "Glu C(bicarbonate)" => [r"E[^E, P]"],          # EX, X≠P/E
    "Glu C(phosphate)" => [r"E[^E, P]", r"D[^E]"],  # EX, X≠P/E, DX, X≠E
    "Asp N/Glu C(bicarbonate)" => [r".D", r"E."],   # Asp N/Glu C(bicarbonate) + EP + EE
    "Lys C/Glu C(phosphate)" => [r"K.", r"E[^E, P]", r"D[^E]"],
    "Chymotrypsin" => [r"[F, Y, W][^P]"],                           # F/Y/WX, X≠P
    "Chymotrypsin(low specificity)" => [r"[F, Y, W, L, X][^P]"],    # F/Y/W/M/LX, X≠P
    "Trypsin/Chymotrypsin" => [r"[F, Y, W, L, X, K, R][^P]"],
    "MAFA" => [r"D."]  # DX
)

const ADDUCT_FN = (
    Dict(
        "[M]" => identity,
        "[M+H]+" => (m, ) -> (m + first(ADD_MS)["H"]), 
        "[M+2H]2+" => (m, ) -> (m / 2 + first(ADD_MS)["H"]),
        "[M-H]-" => (m, ) -> (m - first(ADD_MS)["H"]), 
        "[M-2H]2-" => (m, ) -> (m / 2 - first(ADD_MS)["H"])
    ),
    Dict(
        "[M]" => identity,
        "[M+H]+" => (m, ) -> (m + last(ADD_MS)["H"]), 
        "[M+2H]2+" => (m, ) -> (m / 2 + last(ADD_MS)["H"]),
        "[M-H]-" => (m, ) -> (m - last(ADD_MS)["H"]), 
        "[M-2H]2-" => (m, ) -> (m / 2 - last(ADD_MS)["H"])
    )
)

const UPPER_INDEX = Dict(
    "+" => "$(Char(8314))",
    "-" => "$(Char(8315))",
    "2+" => Char(178) * Char(8314),
    "2-" => Char(178) * Char(8315)
)