
# TODO: 
# - tidy up this file


Base.zero(lv::LorentzVector{CS,T}) where {CS,T} = LorentzVector{CS}(ntuple(i->zero(T),4)...)


Base.:(==)(lv1::LorentzVector, lv2::LorentzVector) = false

function Base.:(==)(lv1::LorentzVector{CS}, lv2::LorentzVector{CS}) where {CS<:LorentzVectorBase.AbstractCoordinateSystem}
    return lv1.data==lv2.data
end

Base.isapprox(lv1::LorentzVector,lv2::LorentzVector; args...) = false

function Base.isapprox(
    lv1::LorentzVector{CS}, 
    lv2::LorentzVector{CS}; 
    atol::Real=0, 
    rtol::Real=atol>0 ? 0 : sqrt(eps())
) where {CS<:LorentzVectorBase.AbstractCoordinateSystem}
    isapprox(lv1.data,lv2.data;atol,rtol)
end


# phi is in [0,2pi) by definition in LorentzVectorBase
function phi02pi(lv::LorentzVector)
    return phi(lv)
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

function fast_mass(v1::LorentzVectorCyl, v2::LorentzVectorCyl)
    # Calculate mass directly. Same as (v1+v2).mass except
    # this skips the intermediate pt, eta, phi calculations.
    # ~4x faster than (v1+v2).mass
    pt1, pt2 = v1.pt, v2.pt
    eta1, eta2 = v1.eta, v2.eta
    phi1, phi2 = v1.phi, v2.phi
    m1, m2 = v1.mass, v2.mass
    
    # note, massless approximation is
    # mass = sqrt(max(2*pt1*pt2*(cosh(eta1-eta2) - cos(phi1-phi2)), zero(pt1)))
    
    sinheta1 = sinh(eta1)
    sinheta2 = sinh(eta2)
    tpt12 = 2*pt1*pt2
    return @fastmath sqrt(max(fma(m1, m1, m2^2)
        + 2*sqrt((pt1^2*(1+sinheta1^2) + m1^2)*(pt2^2*(1+sinheta2^2) + m2^2))
        - tpt12*sinheta1*sinheta2
        - tpt12*cos(phi1-phi2), zero(pt1)))
end
