module LorentzVectorHEPStaticArraysExt

using LorentzVectorHEP, StaticArrays


function StaticArrays.SArray{I,T}(v::LorentzVector) where {I<:Tuple{4},T}
    @SVector T[v.t, v.x, v.y, v.z]
end
function StaticArrays.SArray{I,T}(v::LorentzVectorCyl) where {I<:Tuple{4},T}
    @SVector T[v.pt, v.eta, v.phi, v.mass]
end

#====================================================================================================
NOTE: confusingly, we *are* allowed to set SArray(v) but *NOT* SVector(v).
This is because the latter would be a constructor for a singleton array [v]
====================================================================================================#

for LV ∈ (:LorentzVector, :LorentzVectorCyl)
    @eval StaticArrays.SArray{I}(v::$LV{T}) where {I<:Tuple{4},T} = SArray{I,T}(v)
    
    @eval StaticArrays.SVector{4,T}(v::$LV) where {T} = SArray{Tuple{4},T}(v)
    @eval StaticArrays.SVector{4}(v::$LV{T}) where {T} = SVector{4,T}(v)

    @eval StaticArrays.SArray(v::$LV) = SVector{4}(v)

    @eval function StaticArrays.StaticArray{I,T}(v::$LV) where {I<:Tuple{4},T}
        SArray{Tuple{4},T}(v)
    end
    @eval function StaticArrays.StaticArray{I}(v::$LV) where {I<:Tuple{4}}
        SArray{Tuple{4}}(v)
    end
    @eval StaticArrays.StaticArray(v::$LV) = SArray(v)

    @eval function LorentzVectorHEP.$LV{T}(v::StaticVector{4}) where {T}
        $LV(ntuple(j -> convert(T, v[j]), Val(4))...)
    end
    @eval LorentzVectorHEP.$LV(v::StaticVector{4}) = $LV{eltype(v)}(v)
end


end
