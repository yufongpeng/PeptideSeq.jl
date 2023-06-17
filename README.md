# PeptideSeq

[![Build Status](https://github.com/yufongpeng/PeptideSeq.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yufongpeng/PeptideSeq.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yufongpeng/PeptideSeq.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yufongpeng/PeptideSeq.jl)

*PeptideSeq.jl* is a julia package for predicting peptide sequence after digestion, adding modification and generating expected fragments in mass spectrometry. Digestion enzyme and modification can be customized.

## Installation
This package is not registered yet. Insall it through github:
```julia
julia> using Pkg; Pkg.add("https://github.com/yufongpeng/PeptideSeq.jl")
```

## Example
```julia
julia> using PeptideSeq

julia> p = Protein("DPCHKPKRRKP")
Protein:      DPCHKPKRRKP
Modification: 
Enzyme:       


julia> digest!(p, 2, "Trypsin") # Digest protein with Trypsin
Protein:      DPCHKPKRRKP
Modification: 
Enzyme:       Trypsin
Peptides:
┌──────────────────┬──────────┬─────────┬────────┬───────────────┬──────────────┐
│ Peptide_sequence │ Position │ Mass    │ Adduct │ Miss_cleavage │ Modification │
├──────────────────┼──────────┼─────────┼────────┼───────────────┼──────────────┤
│          DPCHKPK │      1:7 │ 823.401 │ [M]    │             0 │              │
│                R │      8:8 │ 174.112 │ [M]    │             0 │              │
│                R │      9:9 │ 174.112 │ [M]    │             0 │              │
│               KP │    10:11 │ 243.158 │ [M]    │             0 │              │
│         DPCHKPKR │      1:8 │ 979.502 │ [M]    │             1 │              │
│               RR │      8:9 │ 330.213 │ [M]    │             1 │              │
│              RKP │     9:11 │ 399.259 │ [M]    │             1 │              │
│        DPCHKPKRR │      1:9 │  1135.6 │ [M]    │             2 │              │
│             RRKP │     8:11 │ 555.361 │ [M]    │             2 │              │
└──────────────────┴──────────┴─────────┴────────┴───────────────┴──────────────┘

julia> modify!(p, "3NPH") # Derivatize peptides with 3NPH. It can be also modified first and then digested.
Protein:      DPCHKPKRRKP
Modification: 3NPH
Enzyme:       Trypsin
Peptides:
┌──────────────────┬──────────┬─────────┬────────┬───────────────┬───────────────┐
│ Peptide_sequence │ Position │ Mass    │ Adduct │ Miss_cleavage │ Modification  │
├──────────────────┼──────────┼─────────┼────────┼───────────────┼───────────────┤
│          DPCHKPK │      1:7 │ 1093.49 │ [M]    │             0 │ 3NPH (@ 1, 7) │
│                R │      8:8 │ 309.155 │ [M]    │             0 │ 3NPH (@ 8)    │
│                R │      9:9 │ 309.155 │ [M]    │             0 │ 3NPH (@ 9)    │
│               KP │    10:11 │ 378.202 │ [M]    │             0 │ 3NPH (@ 11)   │
│         DPCHKPKR │      1:8 │ 1249.59 │ [M]    │             1 │ 3NPH (@ 1, 8) │
│               RR │      8:9 │ 465.256 │ [M]    │             1 │ 3NPH (@ 9)    │
│              RKP │     9:11 │ 534.303 │ [M]    │             1 │ 3NPH (@ 11)   │
│        DPCHKPKRR │      1:9 │ 1405.69 │ [M]    │             2 │ 3NPH (@ 1, 9) │
│             RRKP │     8:11 │ 690.404 │ [M]    │             2 │ 3NPH (@ 11)   │
└──────────────────┴──────────┴─────────┴────────┴───────────────┴───────────────┘

julia> ionize!(p, "[M+H]+", "[M+2H]2+")
Protein:      DPCHKPKRRKP
Modification: 3NPH       
Enzyme:       Trypsin    
Peptides:
┌──────────────────┬──────────┬─────────┬──────────┬───────────────┬───────────────┐
│ Peptide_sequence │ Position │ Mass    │ Adduct   │ Miss_cleavage │ Modification  │
├──────────────────┼──────────┼─────────┼──────────┼───────────────┼───────────────┤
│          DPCHKPK │      1:7 │  1094.5 │ [M+H]⁺   │             0 │ 3NPH (@ 1, 7) │
│          DPCHKPK │      1:7 │ 547.752 │ [M+2H]²⁺ │             0 │ 3NPH (@ 1, 7) │
│                R │      8:8 │ 310.163 │ [M+H]⁺   │             0 │ 3NPH (@ 8)    │
│                R │      8:8 │ 155.585 │ [M+2H]²⁺ │             0 │ 3NPH (@ 8)    │
│                R │      9:9 │ 310.163 │ [M+H]⁺   │             0 │ 3NPH (@ 9)    │
│                R │      9:9 │ 155.585 │ [M+2H]²⁺ │             0 │ 3NPH (@ 9)    │
│               KP │    10:11 │ 379.209 │ [M+H]⁺   │             0 │ 3NPH (@ 11)   │
│               KP │    10:11 │ 190.109 │ [M+2H]²⁺ │             0 │ 3NPH (@ 11)   │
│         DPCHKPKR │      1:8 │  1250.6 │ [M+H]⁺   │             1 │ 3NPH (@ 1, 8) │
│         DPCHKPKR │      1:8 │ 625.802 │ [M+2H]²⁺ │             1 │ 3NPH (@ 1, 8) │
│               RR │      8:9 │ 466.264 │ [M+H]⁺   │             1 │ 3NPH (@ 9)    │
│               RR │      8:9 │ 233.636 │ [M+2H]²⁺ │             1 │ 3NPH (@ 9)    │
│              RKP │     9:11 │  535.31 │ [M+H]⁺   │             1 │ 3NPH (@ 11)   │
│              RKP │     9:11 │ 268.159 │ [M+2H]²⁺ │             1 │ 3NPH (@ 11)   │
│        DPCHKPKRR │      1:9 │  1406.7 │ [M+H]⁺   │             2 │ 3NPH (@ 1, 9) │
│        DPCHKPKRR │      1:9 │ 703.853 │ [M+2H]²⁺ │             2 │ 3NPH (@ 1, 9) │
│             RRKP │     8:11 │ 691.412 │ [M+H]⁺   │             2 │ 3NPH (@ 11)   │
│             RRKP │     8:11 │  346.21 │ [M+2H]²⁺ │             2 │ 3NPH (@ 11)   │
└──────────────────┴──────────┴─────────┴──────────┴───────────────┴───────────────┘

julia> fragmentation(p.peptides[16]) # Fragmentation of 16th peptide
Peptide:      DPCHKPKRR
Modification: 3NPH (@ 1, 9)
Precursor:    703.8527195000001
Adduct:       [M+2H]²⁺
Fragment ions:    
┌──────┬─────────┐
│ Type │ Mass    │
├──────┼─────────┤
│ b₁⁺  │ 251.078 │
│ b₂⁺  │ 348.131 │
│ b₃⁺  │  451.14 │
│ b₄⁺  │ 588.199 │
│ b₅⁺  │ 716.294 │
│ b₆⁺  │ 813.347 │
│ b₆²⁺ │ 407.177 │
│ b₇⁺  │ 941.442 │
│ b₇²⁺ │ 471.225 │
│ b₈⁺  │ 1097.54 │
│ b₈²⁺ │ 549.275 │
│ y₁⁺  │ 310.163 │
│ y₂⁺  │ 466.264 │
│ y₃⁺  │ 594.359 │
│ y₄⁺  │ 691.412 │
│ y₅⁺  │ 819.507 │
│ y₆⁺  │ 956.565 │
│ y₆²⁺ │ 478.787 │
│ y₇⁺  │ 1059.57 │
│ y₇²⁺ │ 530.291 │
│ y₈⁺  │ 1156.63 │
│ y₈²⁺ │ 578.818 │
└──────┴─────────┘

julia> data = read_msdial("test/data.txt") # MSDIAL data
Table with 6 columns and 411 rows:
      Alignment ID  Average Rt(min)  Average Mz  Metabolite name  Adduct type  MS/MS spectrum
    ┌────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 1  │ 0             0.691            100.077     Unknown          [M+H-2H2O]+  [50.1298 0.0; 50.3013 0.0; 50.3086 0.0; 50.3209 0.0; 50.3…
 2  │ 1             15.696           102.128     Unknown          [M+H]+       [50.3038 0.0; 50.5715 1.0; 50.9706 1.0; 51.0595 0.0; 51.1…
 3  │ 2             19.7             102.129     Unknown          [M+H]+       [50.147 0.0; 50.2033 0.0; 50.2523 0.0; 50.3038 0.0; 50.38…
 4  │ 3             16.241           102.129     Unknown          [M+H]+       Float64[]
 5  │ 4             17.278           102.129     Unknown          [M+H]+       [50.1228 0.0; 50.2895 0.0; 50.4195 0.0; 50.5325 0.0; 50.5…
 6  │ 5             4.481            118.085     Unknown          [M+H]+       [50.3526 0.0; 50.4704 0.0; 50.5515 0.0; 50.5663 0.0; 50.6…
 7  │ 6             11.167           118.086     Unknown          [M+H]+       Float64[]
 8  │ 7             0.682            132.103     NAE 4:0          [M+H]+       [50.064 0.0; 50.1986 0.0; 50.2673 0.0; 50.2746 0.0; 50.28…
 9  │ 8             3.383            149.026     Unknown          [M+H+2Na]+   [50.0222 0.0; 50.1054 0.0; 50.2719 0.0; 50.3332 0.0; 50.4…
 10 │ 9             6.261            149.026     Unknown          [M+H-H2O]+   [50.154 0.0; 50.3697 0.0; 50.4311 0.0; 50.4925 0.0; 50.68…
 11 │ 10            0.605            152.02      Unknown          [M+H]+       [50.227 0.0; 50.5264 0.0; 50.5363 0.0; 50.7774 0.0; 50.98…
 ⋮  │      ⋮               ⋮             ⋮              ⋮              ⋮                                   ⋮


julia> result = find_peptides(data, db) # Find peptides in data
Table with 11 columns and 411 rows:
      Alignment ID  Average Rt(min)  Average Mz  Metabolite name  Adduct type  MS/MS spectrum        Peptide Matched  Peptide  ⋯
    ┌───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 1  │ 0             0.691            100.077     Unknown          [M+H-2H2O]+  [50.1298 0.0; 50.30…  false            nothing  ⋯
 2  │ 1             15.696           102.128     Unknown          [M+H]+       [50.3038 0.0; 50.57…  false            nothing  ⋯
 3  │ 2             19.7             102.129     Unknown          [M+H]+       [50.147 0.0; 50.203…  false            nothing  ⋯
 4  │ 3             16.241           102.129     Unknown          [M+H]+       Float64[]             false            nothing  ⋯
 5  │ 4             17.278           102.129     Unknown          [M+H]+       [50.1228 0.0; 50.28…  false            nothing  ⋯
 6  │ 5             4.481            118.085     Unknown          [M+H]+       [50.3526 0.0; 50.47…  false            nothing  ⋯
 7  │ 6             11.167           118.086     Unknown          [M+H]+       Float64[]             false            nothing  ⋯
 8  │ 7             0.682            132.103     NAE 4:0          [M+H]+       [50.064 0.0; 50.198…  false            nothing  ⋯
 9  │ 8             3.383            149.026     Unknown          [M+H+2Na]+   [50.0222 0.0; 50.10…  false            nothing  ⋯
 10 │ 9             6.261            149.026     Unknown          [M+H-H2O]+   [50.154 0.0; 50.369…  false            nothing  ⋯
 11 │ 10            0.605            152.02      Unknown          [M+H]+       [50.227 0.0; 50.526…  false            nothing  ⋯
 ⋮  │      ⋮               ⋮             ⋮              ⋮              ⋮                ⋮                   ⋮            ⋮     ⋱

julia> result.Peptide[121]
Peptide:      DPCHKPKR
Modification: 3NPH (@ 1)
Precursor:    558.2805335
Adduct:       [M+2H]²⁺
Fragment ions:
┌──────┬─────────┐
│ Type │ Mass    │
├──────┼─────────┤
│ b₁⁺  │ 251.078 │
│ b₂⁺  │ 348.131 │
│ b₃⁺  │  451.14 │
│ b₄⁺  │ 588.199 │
│ b₅⁺  │ 716.294 │
│ b₆⁺  │ 813.347 │
│ b₆²⁺ │ 407.177 │
│ b₇⁺  │ 941.442 │
│ b₇²⁺ │ 471.225 │
│ y₁⁺  │ 175.119 │
│ y₂⁺  │ 303.214 │
│ y₃⁺  │ 400.267 │
│ y₄⁺  │ 528.362 │
│ y₅⁺  │ 665.421 │
│ y₆⁺  │  768.43 │
│ y₆²⁺ │ 384.719 │
│ y₇⁺  │ 865.483 │
│ y₇²⁺ │ 433.245 │
└──────┴─────────┘

```