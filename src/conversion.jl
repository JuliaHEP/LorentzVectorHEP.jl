
# TODO: 
# - exploit the Base.convert instead of calling constructors
# - unify this implementation for every coordinate system by using coordinate_names


# conversion
function LorentzVector{CS}(lv::LorentzVector) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}
    return LorentzVector(
        LorentzVectorBase.px(lv), 
        LorentzVectorBase.py(lv), 
        LorentzVectorBase.pz(lv),
        LorentzVectorBase.pt(lv), 
    )
end

function LorentzVectorCyl(lv::LorentzVector)
    return LorentzVectorCyl(
        LorentzVectorBase.pt(lv),
        LorentzVectorBase.eta(lv),
        LorentzVectorBase.phi(lv),
        LorentzVectorBase.mass(lv),
    )
end

function fromPtEtaPhiE(pt, eta, phi, E) 
    m2 = E^2 - pt^2 - (sinh(eta) * pt)^2
    m = sign(m2) * sqrt(abs(m2)) 
    return LorentzVectorCyl(pt, eta, phi, m) 
end

function fromPxPyPzM(px, py, pz, m)
    e = sqrt(muladd(px, px, py^2) + muladd(pz, pz, m^2))
    return LorentzVector(e, px, py, pz)
end
