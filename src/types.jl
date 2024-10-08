
# TODO: 
# - get AbstractLorentzVector from LorentzVectors.jl
#   after https://github.com/JLTastet/LorentzVectors.jl/pull/12 is merged
abstract type AbstractLorentzVector end

struct LorentzVector{CS<:LorentzVectorBase.AbstractCoordinateSystem,T} <:AbstractLorentzVector
    coord_sys::CS
    data::NTuple{4,T}

    LorentzVector{CS}(comp1::T,comp2::T,comp3::T,comp0::T) where {CS,T} = new{CS,T}(CS(),(comp1,comp2,comp3,comp0))

    function LorentzVector(cs::CS,comp1::T,comp2::T,comp3::T,comp0::T) where {CS<:LorentzVectorBase.AbstractCoordinateSystem,T}

        new{CS,T}(cs,(comp1,comp2,comp3,comp0))
    end
end

# type promotion for components
LorentzVector(cs::LorentzVectorBase.AbstractCoordinateSystem,c1,c2,c3,c0) = LorentzVector(cs,promote(c1,c2,c3,c0)...)
LorentzVector{CS}(c1,c2,c3,c0) where CS = LorentzVector{CS}(promote(c1,c2,c3,c0)...) 

# default coordinate system
LorentzVector(c1,c2,c3,c0) = LorentzVector(LorentzVectorBase.XYZT(),c1,c2,c3,c0)

Base.broadcastable(lv::LorentzVector) = Ref(lv)

# all supported coordinate system types
const SUPPORTED_COORDINATE_SYSTEMS = InteractiveUtils.subtypes(LorentzVectorBase.AbstractCoordinateSystem)

# implementation of the kinematic interface
LorentzVectorBase.coordinate_system(lv::LorentzVector{CS}) where {CS<:LorentzVectorBase.AbstractCoordinateSystem}= lv.coord_sys

# For each possible coordiante system, this implements the respective getter functions.
# The order of the components is given by coordinate_names(cs) for every coodinate system
# cs. 
for cs in SUPPORTED_COORDINATE_SYSTEMS
    for (i,func) in enumerate(coordinate_names(cs()))
        eval(
            quote
                (LorentzVectorBase.$func)(lv::LorentzVector{$cs}) = lv.data[$i]
            end,
        )
    end
end

# type alias for cylindrical coordinates
const LorentzVectorCyl = LorentzVector{LorentzVectorBase.PtEtaPhiM}

