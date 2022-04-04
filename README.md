# LorentzVectorHEP

```julia
pt(lv::LorentzVectorCyl) = lv.pt
eta(lv::LorentzVectorCyl) = lv.eta
phi(lv::LorentzVectorCyl) = lv.phi
mass(lv::LorentzVectorCyl) = lv.mass
px(v::LorentzVectorCyl) = v.pt * cos(v.phi)
py(v::LorentzVectorCyl) = v.pt * sin(v.phi)
pz(v::LorentzVectorCyl) = v.pt * sinh(v.eta)
energy(v::LorentzVectorCyl) = sqrt(px(v)^2 + py(v)^2 + pz(v)^2 + v.mass^2)
```
