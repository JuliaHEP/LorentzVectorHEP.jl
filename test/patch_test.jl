
# These patches restore compatibility between the new implementation with the old tests. 
# This can be removed, if the tests are refactored as well. 

# adopting the new order of arguments from LorentzVectorBase
old_LorentzVector(t,x,y,z) = LorentzVector(x,y,z,t) 

# adopting the new conversion syntax
old_LorentzVector(lv::LorentzVectorCyl) = LorentzVector{LorentzVectorBase.XYZT}(lv) 

# overwrite the getproperty for the elementary accessors
function Base.getproperty(lv::LorentzVector,sym::Symbol)
    if sym in (:x,:y,:z,:t,:phi,:pt,:eta,:mass)
        return eval(
            quote
                LorentzVectorBase.$sym($lv)
            end
        )
    else
        return getfield(lv,sym)
    end
end
