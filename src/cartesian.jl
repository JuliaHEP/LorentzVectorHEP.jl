Base.zero(lv::T) where T<:LorentzVector = T(0,0,0,0)
mass2(lv::LorentzVector) = dot(lv, lv)
mass(lv::LorentzVector) = mass2(lv) < 0.0 ? sqrt(-mass2(lv)) : sqrt(mass2(lv))
pt2(lv::LorentzVector) = muladd(lv.x, lv.x, lv.y^2)
pt(lv::LorentzVector) = sqrt(pt2(lv))
mt2(lv::LorentzVector) = lv.t^2 - lv.z^2
mt(lv::LorentzVector) = mt2(lv)<0 ? -sqrt(-mt2(lv)) : sqrt(mt2(lv))
mag(lv::LorentzVector) = sqrt(muladd(lv.x, lv.x, lv.y^2) + lv.z^2)
energy(lv::LorentzVector) = lv.t
px(lv::LorentzVector) = lv.x
py(lv::LorentzVector) = lv.y
pz(lv::LorentzVector) = lv.z

@inline function CosTheta(lv::LorentzVector)
    fZ = lv.z
    ptot = mag(lv)
    return ifelse(ptot == 0.0, 1.0, fZ / ptot)
end

"""Pseudorapidity"""
function eta(lv::LorentzVector)
    cosTheta = CosTheta(lv)
    (cosTheta^2 < 1.0) && return -0.5 * log((1.0 - cosTheta) / (1.0 + cosTheta))
    fZ = lv.z
    iszero(fZ) && return 0.0
    # Warning("PseudoRapidity","transverse momentum = 0! return +/- 10e10");
    fZ > 0.0 && return 10e10
    return -10e10
end
const η = eta

"""Rapidity"""
function rap(lv::LorentzVector)
    pt_squared = pt2(lv)
    abspz = abs(pz(lv))
    if (energy(lv) == abspz) && (pt_squared == 0.0)
        return (-1)^(pz(lv) < 0)*(1e5 + abspz) # a very large value that depends on pz
    end
    m2 = max((energy(lv) + pz(lv))*(energy(lv) - pz(lv)) - pt_squared, 0.0) # mass^2
    E_plus_z = energy(lv) + abspz
    return (-1)^(pz(lv) > 0) * 0.5*log((pt_squared + m2)/(E_plus_z^2))
end
# Don't export "y" as an alias as for a normal 4-vector it's the second space coordinate

function phi(lv::LorentzVector)
    return (lv.x == 0.0 && lv.y == 0.0) ? 0.0 : atan(lv.y, lv.x)
end
const ϕ = phi

function phi02pi(lv::LorentzVector)
    return phi(lv) < 0.0 ? phi(lv) + 2π : phi(lv)
end

function phi_mpi_pi(x)
    twopi = 2pi
    while (x >= pi)
        x -= twopi
    end
    while (x < -pi)
        x += twopi
    end
    return x
end

deltaeta(lv1, lv2) = eta(lv1) - eta(lv2)

deltaphi(lv1, lv2) = phi_mpi_pi(phi(lv1) - phi(lv2))

function deltar(lv1, lv2)
    deta = eta(lv1) - eta(lv2)
    dphi = deltaphi(lv1, lv2)
    return sqrt(fma(deta, deta, dphi * dphi))
end

function fromPxPyPzM(px, py, pz, m)
    e = sqrt(muladd(px, px, py^2) + muladd(pz, pz, m^2))
    return LorentzVector(e, px, py, pz)
end
