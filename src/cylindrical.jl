struct LorentzVectorCyl{T <: AbstractFloat}
    pt::T
    eta::T
    phi::T
    mass::T
end

Base.show(io::IO, v::LorentzVectorCyl) = print(io, "$(typeof(v))(pt=$(v.pt), eta=$(v.eta), phi=$(v.phi), mass=$(v.mass))")

Base.broadcastable(v::LorentzVectorCyl) = Ref(v)

"""
    zero(LorentzVectorCyl)
    zero(v::LorentzVectorCyl{T})
Constructs a zero four-vector.
"""
Base.zero(::Type{LorentzVectorCyl{T}}) where T = LorentzVectorCyl{T}(zero(T), zero(T), zero(T), zero(T)) 
Base.zero(::Type{LorentzVectorCyl}) = zero(LorentzVectorCyl{Float64})
Base.zero(lv::T) where T<:LorentzVectorCyl = T(0,0,0,0)

pt(lv::LorentzVectorCyl) = lv.pt
eta(lv::LorentzVectorCyl) = lv.eta
phi(lv::LorentzVectorCyl) = lv.phi
mass(lv::LorentzVectorCyl) = lv.mass
px(v::LorentzVectorCyl) = v.pt * cos(v.phi)
py(v::LorentzVectorCyl) = v.pt * sin(v.phi)
pz(v::LorentzVectorCyl) = v.pt * sinh(v.eta)
energy(v::LorentzVectorCyl) = sqrt(px(v)^2 + py(v)^2 + pz(v)^2 + v.mass^2)

function Base.:*(v::LorentzVectorCyl{T}, k::Real) where T
    LorentzVectorCyl{T}(v.pt*k, v.eta, v.phi, v.mass*k)
end

function Base.:+(v1::LorentzVectorCyl{T}, v2::LorentzVectorCyl{W}) where {T,W}
    m1, m2 = max(v1.mass, zero(v1.pt)), max(v2.mass, zero(v2.pt))
    
    px1, px2 = px(v1), px(v2)
    py1, py2 = py(v1), py(v2)
    pz1, pz2 = pz(v1), pz(v2)
    e1 = sqrt(px1^2 + py1^2 + pz1^2 + m1^2)
    e2 = sqrt(px2^2 + py2^2 + pz2^2 + m2^2)
    
    sumpx = px1+px2
    sumpy = py1+py2
    sumpz = pz1+pz2
    
    ptsq = sumpx^2 + sumpy^2
    pt = sqrt(ptsq)
    eta = asinh(sumpz/pt)
    phi = atan(sumpy, sumpx)
    mass = sqrt(max(muladd(m1, m1, m2^2) + 2*e1*e2 - 2*(muladd(px1, px2, py1*py2) + pz1*pz2), zero(v1.pt)))
    return LorentzVectorCyl(pt,eta,phi,mass)
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


# https://root.cern.ch/doc/v606/GenVector_2VectorUtil_8h_source.html#l00061
@inline function deltaphi(v1::LorentzVectorCyl, v2::LorentzVectorCyl)
    dphi = v2.phi - v1.phi
    dphi = dphi - 2*pi*(dphi > pi) + 2*pi*(dphi <= -pi)
    return dphi
end
const Δϕ = deltaphi

@inline function deltaeta(v1::LorentzVectorCyl, v2::LorentzVectorCyl)
    return v2.eta - v1.eta
end
const Δη = deltaeta

@inline function deltar2(v1::LorentzVectorCyl, v2::LorentzVectorCyl)
    dphi = deltaphi(v1,v2)
    deta = deltaeta(v1,v2)
    return muladd(dphi, dphi, deta^2)
end
deltar(v1::LorentzVectorCyl, v2::LorentzVectorCyl) = sqrt(deltar2(v1, v2))
const ΔR = deltar


function fromPtEtaPhiE(pt, eta, phi, E) 
    m2 = E^2 - pt^2 - (sinh(eta) * pt)^2
    m = sign(m2) * sqrt(abs(m2)) 
    return LorentzVectorCyl(pt, eta, phi, m) 
end
  
