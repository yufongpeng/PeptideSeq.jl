# PeptideSeq

[![Build Status](https://github.com/yufongpeng/PeptideSeq.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yufongpeng/PeptideSeq.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yufongpeng/PeptideSeq.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yufongpeng/PeptideSeq.jl)

*PeptideSeq.jl* is a julia package for predicting digested peptide sequence, adding modification and generating expected fragments in mass spectrometry. 

## Installation
This package is not yet registered. Insall it through github:
```julia
julia> using Pkg; Pkg.add("https://github.com/yufongpeng/PeptideSeq.jl")
```

## Example
```julia
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
```