module LorentzVectorHEP

using LinearAlgebra

export LorentzVectorCyl, LorentzVector


export px, py, pz, energy, fast_mass, pt, rapidity, eta, phi, theta, mass
export deltaphi, deltar, deltaeta, deltatheta
export ΔR, Δϕ, Δη
export fromPtEtaPhiE


"""
    LorentzVector(t, x, y, z)

Lorentz 4-vector, as used in Special Relativity.

The metric convention is g = diag(+1,-1,-1,-1). No distinction is made between
co- and contravariant vectors.
"""
struct LorentzVector{T <: AbstractFloat}
    t :: T
    x :: T
    y :: T
    z :: T
end

"""
    LorentzVector(t, x, y, z)

Promoting constructors for LorentzVector{T}.
"""
LorentzVector(t, x, y, z) = LorentzVector(promote(t, x, y, z)...)
LorentzVector(t::T, x::T, y::T, z::T) where {T <: Union{Integer, Rational, Irrational}} =
    LorentzVector(float(t), x, y, z)

"""
    dot(u, v)
    u⋅v

Inner product of 4-vectors, in the Minkowsky metric (+,-,-,-).
"""
function LinearAlgebra.dot(u::LorentzVector, v::LorentzVector)
    @fastmath u.t*v.t - u.x*v.x - u.y*v.y - u.z*v.z
end


function Base.:(+)(u::LorentzVector, v::LorentzVector)
    @fastmath LorentzVector(u.t + v.t, u.x + v.x, u.y + v.y, u.z + v.z)
end

function Base.:(-)(u::LorentzVector)
    @fastmath LorentzVector(-u.t, -u.x, -u.y, -u.z)
end

function Base.:(-)(u::LorentzVector, v::LorentzVector)
    @fastmath u + (-v)
end

function Base.:(*)(λ::Number, u::LorentzVector)
    @fastmath LorentzVector(λ*u.t, λ*u.x, λ*u.y, λ*u.z)
end

function Base.:(*)(u::LorentzVector, λ::Number)
    @fastmath λ * u
end

function Base.:(/)(u::LorentzVector, λ::Number)
    @fastmath u * (one(λ) / λ)
end

function Base.:(==)(u::LorentzVector, v::LorentzVector)
    u.t == v.t && u.x == v.x && u.y == v.y && u.z == v.z
end

function Base.isapprox(u::LorentzVector, v::LorentzVector;
                  atol::Real=0,
                  rtol::Real=atol>0 ? 0 : √min(eps(typeof(u.x)), eps(typeof(v.x))))

    err = max(abs(u.t - v.t), sqrt((u.x - v.x)^2 + (u.y - v.y)^2 + (u.z - v.z)^2))
    err <= max(atol, rtol*max(abs(u.t), abs(v.t), sqrt(u.x^2 + u.y^2 + u.z^2), sqrt(v.x^2 + v.y^2 + v.z^2)))
end

include("cartesian.jl")
include("cylindrical.jl")

const Δϕ = deltaphi
const Δη = deltaeta
const ΔR = deltar
const Δθ = deltatheta

export ΔR, Δϕ, Δη, Δθ

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
