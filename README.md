# LorentzVectorHEP

Provides two types (and the conversion between the two):

- `LorentzVector` (energy, px, py, pz) (from [LorentzVectors.jl](https://github.com/JLTastet/LorentzVectors.jl))
- `LorentzVectorCyl` (pt, eta, phi, mass)

you can also use `fromPtEtaPhiE(pt, eta, phi, energy) --> LorentzVectorCyl`.

and these functions for both of them:

```julia
px, py, pz, energy, pt, rapidity, eta, phi, mass
```

as well as these utility functions:

```julia
deltar, deltaphi, deltaeta, mt, mt2
```

(some of them have aliases, `ΔR, Δϕ, Δη`)

There are some unexported methods which are useful for more specialist use cases:

```julia
mass2, pt2, mt, mt2, mag
```

## LHC coordinate system

![LHC Coordinate System](https://cds.cern.ch/record/1699952/files/Figures_T_Coordinate.png)
