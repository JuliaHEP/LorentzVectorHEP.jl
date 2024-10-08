function Base.show(io::IO, ::MIME"text/plain", lv::LorentzVector{CS}) where CS
    println(io, "LorentzVector: $(nameof(CS))")
    for comp in coordinate_names(lv.coord_sys)
        println(io, "        $comp = $(eval(
			quote
				LorentzVectorBase.$comp($lv)
			end
			))")
    end
    return nothing
end

function Base.show(io::IO,  lv::LorentzVector{CS}) where CS
	cs_name = nameof(CS)
	ret_str = "LorentzVector{$cs_name}("
	for comp in coordinate_names(lv.coord_sys)
        ret_str *= "$comp = $(eval(
			quote
				LorentzVectorBase.$comp($lv)
			end
			)), "
    end
	ret_str = ret_str[1:end-2]*")"

	println(ret_str)
	return nothing
end
