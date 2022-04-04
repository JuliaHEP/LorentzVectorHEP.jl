# LorentzVectorHEP

Provides two types (and the conversion between the two):
- `LorentzVector` (energy, px, py, pz) (from [LorentzVectors.jl](https://github.com/JLTastet/LorentzVectors.jl))
- `LorentzVectorCyl` (pt, eta, phi, mass)

and these functions for both of them:
```julia
px, py, pz, energy, pt, eta, phi, mass
```


as well as these utility functions:
```julia
deltar, deltaeta, mt, mt2
```


## LHC coordinate system

![](https://cds.cern.ch/record/1699952/files/Figures_T_Coordinate.png)
