# TODO: 
# - think about implementing arithmetic for all possbile coordinate systems by converting
# to cartesian, perform the operation, and then converting back to the original
# coordiantes.


### cartesian coordinates

function Base.:(+)(lv1::LorentzVector{CS},lv2::LorentzVector{CS}) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}

    return LorentzVector{CS}(
        LorentzVectorBase.x(lv1)+LorentzVectorBase.x(lv2),
        LorentzVectorBase.y(lv1)+LorentzVectorBase.y(lv2),
        LorentzVectorBase.z(lv1)+LorentzVectorBase.z(lv2),
        LorentzVectorBase.t(lv1)+LorentzVectorBase.t(lv2),
    )
end

function Base.:(-)(lv1::LorentzVector{CS},lv2::LorentzVector{CS}) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}

    return LorentzVector{CS}(
        LorentzVectorBase.x(lv1)-LorentzVectorBase.x(lv2),
        LorentzVectorBase.y(lv1)-LorentzVectorBase.y(lv2),
        LorentzVectorBase.z(lv1)-LorentzVectorBase.z(lv2),
        LorentzVectorBase.t(lv1)-LorentzVectorBase.t(lv2),
    )
end

function Base.:(-)(lv1::LorentzVector{CS}) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}
    return LorentzVector{CS}(
        -LorentzVectorBase.x(lv1),
        -LorentzVectorBase.y(lv1),
        -LorentzVectorBase.z(lv1),
        -LorentzVectorBase.t(lv1)
    )
end


function Base.:(*)(a::Number,lv1::LorentzVector{CS}) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}

    return LorentzVector{CS}(
        a*LorentzVectorBase.x(lv1),
        a*LorentzVectorBase.y(lv1),
        a*LorentzVectorBase.z(lv1),
        a*LorentzVectorBase.t(lv1),
    )
end

function Base.:(*)(lv1::LorentzVector{CS},a::Number) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}

    return LorentzVector{CS}(
        LorentzVectorBase.x(lv1)*a,
        LorentzVectorBase.y(lv1)*a,
        LorentzVectorBase.z(lv1)*a,
        LorentzVectorBase.t(lv1)*a,
    )
end

function Base.:(/)(lv::LorentzVector{CS}, a::Number) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}
    return LorentzVector{CS}(
        LorentzVectorBase.x(lv)/a,
        LorentzVectorBase.y(lv)/a,
        LorentzVectorBase.z(lv)/a,
        LorentzVectorBase.t(lv)/a,
    )
end

function LinearAlgebra.dot(lv1::LorentzVector{CS},lv2::LorentzVector{CS}) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}
    return LorentzVectorBase.t(lv1)*LorentzVectorBase.t(lv2)
        - LorentzVectorBase.x(lv1)*LorentzVectorBase.x(lv2)
        - LorentzVectorBase.y(lv1)*LorentzVectorBase.y(lv2)
        - LorentzVectorBase.z(lv1)*LorentzVectorBase.z(lv2)
end

### cylindrical coordinates

function Base.:(+)(lv1::LorentzVectorCyl, lv2::LorentzVectorCyl)
    lv1_cart = LorentzVector{LorentzVectorBase.XYZT}(lv1)
    lv2_cart = LorentzVector{LorentzVectorBase.XYZT}(lv2)
    return LorentzVectorCyl(lv1_cart + lv2_cart)
end

function Base.:(-)(lv1::LorentzVectorCyl, lv2::LorentzVectorCyl)
    lv1_cart = LorentzVector{LorentzVectorBase.XYZT}(lv1)
    lv2_cart = LorentzVector{LorentzVectorBase.XYZT}(lv2)
    return LorentzVectorCyl(lv1_cart - lv2_cart)
end

function Base.:(-)(lv::LorentzVectorCyl) 
    lv_cart = LorentzVector{LorentzVectorBase.XYZT}(lv)
    return LorentzVectorCyl(-lv_cart)
end

function Base.:(*)(a::Number,lv::LorentzVectorCyl)
    LorentzVectorCyl(
        a*LorentzVectorBase.pt(lv),
        LorentzVectorBase.eta(lv),
        LorentzVectorBase.phi(lv),
        a*LorentzVectorBase.mass(lv))
end

function Base.:(*)(lv::LorentzVectorCyl, a::Number) 
    LorentzVectorCyl(
        LorentzVectorBase.pt(lv)*a,
        LorentzVectorBase.eta(lv),
        LorentzVectorBase.phi(lv),
        LorentzVectorBase.mass(lv)*a)
end
