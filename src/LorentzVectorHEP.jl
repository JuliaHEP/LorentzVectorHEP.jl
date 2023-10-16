module LorentzVectorHEP

using LorentzVectors # provides x, y, z, t

export LorentzVectorCyl, LorentzVector

export px, py, pz, energy, fast_mass, pt, rapidity, eta, phi, mass
export deltaphi, deltar, deltaeta
export ΔR, Δϕ, Δη
export fromPtEtaPhiE

include("cartesian.jl")
include("cylindrical.jl")

# conversion
function LorentzVector(v::LorentzVectorCyl)
    x = px(v)
    y = py(v)
    z = pz(v)
    t = energy(v)
    return LorentzVector(t, x, y, z)
end

function LorentzVectorCyl(v::LorentzVector)
    t, x, y, z = v.t, v.x, v.y, v.z
    pt2 = muladd(x, x, y^2)
    pt = sqrt(pt2)
    eta = asinh(z/pt)
    phi = atan(y, x)
    mass = sqrt(max(t^2 - pt2 - z^2, 0))
    return LorentzVectorCyl(pt, eta, phi, mass)
end

end
