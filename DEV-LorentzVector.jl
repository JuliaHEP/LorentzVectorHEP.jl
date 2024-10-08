### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ bff0ea06-8574-11ef-0913-f59c3e6f6ee3
begin
	import Pkg
	Pkg.activate(".")


	using LorentzVectorBase
	using BenchmarkTools
end

# ╔═╡ c0d3c40f-6c51-4c29-b1db-5b1e265d8cc8
md"# The type"

# ╔═╡ f1986aa8-5a6f-4643-8bc1-9e65555470e1
begin
	struct LorentzVector{CS<:LorentzVectorBase.AbstractCoordinateSystem,T}
		coord_sys::CS
		data::NTuple{4,T}

		LorentzVector{CS}(comp1::T,comp2::T,comp3::T,comp0::T) where {CS,T} = new{CS,T}(CS(),(comp1,comp2,comp3,comp0))
		
		function LorentzVector(cs::CS,comp1::T,comp2::T,comp3::T,comp0::T) where {CS<:LorentzVectorBase.AbstractCoordinateSystem,T}

			new{CS,T}(cs,(comp1,comp2,comp3,comp0))
		end
	end

	LorentzVector(cs::LorentzVectorBase.AbstractCoordinateSystem,c1,c2,c3,c0) = LorentzVector(cs,promote(c1,c2,c3,c0)...)
	LorentzVector{CS}(c1,c2,c3,c0) where CS = LorentzVector{CS}(promote(c1,c2,c3,c0)...) 
	LorentzVector(c1,c2,c3,c0) = LorentzVector(LorentzVectorBase.XYZT(),c1,c2,c3,c0)

	
	const SUPPORTED_COORDINATE_SYSTEMS = subtypes(LorentzVectorBase.AbstractCoordinateSystem)

	LorentzVectorBase.coordinate_system(lv::LorentzVector{CS}) where {CS<:LorentzVectorBase.AbstractCoordinateSystem}= lv.coord_sys
	
	for cs in SUPPORTED_COORDINATE_SYSTEMS
		for (i,func) in enumerate(coordinate_names(cs()))
			eval(
			    quote
			      (LorentzVectorBase.$func)(lv::LorentzVector{$cs}) = lv.data[$i]
				  end,
			  )
		end
	end


	const LorentzVectorCyl = LorentzVector{LorentzVectorBase.PtEtaPhiM}
	
end

# ╔═╡ 74e9c8c0-2049-4381-b128-22c3326bf588
csys = LorentzVectorBase.PxPyPzE()


# ╔═╡ 9076a01f-63d8-4312-8969-3d874811b82e
lvec = LorentzVector(csys,1.0,0,0,0)

# ╔═╡ 3192a5e1-a59f-4db9-a3c7-2d1fce273a03
coordinate_system(lvec)

# ╔═╡ b83b54dc-e758-4cb7-ae9c-d899d4dbbc05
LorentzVectorBase.mass(lvec)

# ╔═╡ 28934d77-768b-4a37-aa36-c2eb3504ac6d
lvec_cyl = LorentzVectorCyl(0.2,0.1,0,1.0)

# ╔═╡ 2a08d509-bc54-40f1-bf84-b02277f16f70
@benchmark LorentzVectorBase.px($lvec_cyl)

# ╔═╡ c8fe4574-5473-4068-aaaa-3053f82e6899
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

# ╔═╡ 2b2dad33-59ff-4b6c-8522-7ebdc6ac68b7
lvec_cyl

# ╔═╡ 4b349f2b-c197-4e84-a7ff-d7b9114ef1f2
println(lvec,lvec)

# ╔═╡ 3b4822e4-5d27-4a1d-a093-2bb523608b49
md"# Arithmetics"

# ╔═╡ 4c81a5a7-9e34-40b6-b9a5-e8138ad061bf
begin

	function Base.:+(lv1::LorentzVector{CS},lv2::LorentzVector{CS}) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}

		return LorentzVector{CS}(
			LorentzVectorBase.x(lv1)+LorentzVectorBase.x(lv2),
			LorentzVectorBase.y(lv1)+LorentzVectorBase.y(lv2),
			LorentzVectorBase.z(lv1)+LorentzVectorBase.z(lv2),
			LorentzVectorBase.t(lv1)+LorentzVectorBase.t(lv2),
			)
	end

	function Base.:-(lv1::LorentzVector{CS},lv2::LorentzVector{CS}) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}

		return LorentzVector{CS}(
			LorentzVectorBase.x(lv1)-LorentzVectorBase.x(lv2),
			LorentzVectorBase.y(lv1)-LorentzVectorBase.y(lv2),
			LorentzVectorBase.z(lv1)-LorentzVectorBase.z(lv2),
			LorentzVectorBase.t(lv1)-LorentzVectorBase.t(lv2),
			)
	end

	function Base.:*(a::Number,lv1::LorentzVector{CS}) where {CS<:Union{LorentzVectorBase.XYZT,LorentzVectorBase.PxPyPzE}}

		return LorentzVector{CS}(
			a*LorentzVectorBase.x(lv1),
			a*LorentzVectorBase.y(lv1),
			a*LorentzVectorBase.z(lv1),
			a*LorentzVectorBase.t(lv1),
			)
	end
	
	
end

# ╔═╡ aabe344d-c462-4140-bf99-f52086f5cec7
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


# ╔═╡ 662deff8-e425-4442-b524-7c656aeecbd6
3*lvec

# ╔═╡ 54c50afe-314b-46d4-a71d-ded8f028c1e7
lvec-lvec

# ╔═╡ 273b511b-330f-43dd-bf01-99bc980b9fd7
T = Float32

# ╔═╡ f2e150ff-b33f-4697-bebd-790a6437fdf0
n = ntuple(i->zero(T),4)

# ╔═╡ 448fdb1e-7a3c-4105-bd9b-f18d109986a4


# ╔═╡ Cell order:
# ╠═bff0ea06-8574-11ef-0913-f59c3e6f6ee3
# ╠═c0d3c40f-6c51-4c29-b1db-5b1e265d8cc8
# ╠═f1986aa8-5a6f-4643-8bc1-9e65555470e1
# ╠═74e9c8c0-2049-4381-b128-22c3326bf588
# ╠═9076a01f-63d8-4312-8969-3d874811b82e
# ╠═3192a5e1-a59f-4db9-a3c7-2d1fce273a03
# ╠═b83b54dc-e758-4cb7-ae9c-d899d4dbbc05
# ╠═28934d77-768b-4a37-aa36-c2eb3504ac6d
# ╠═2a08d509-bc54-40f1-bf84-b02277f16f70
# ╠═c8fe4574-5473-4068-aaaa-3053f82e6899
# ╠═aabe344d-c462-4140-bf99-f52086f5cec7
# ╠═2b2dad33-59ff-4b6c-8522-7ebdc6ac68b7
# ╠═4b349f2b-c197-4e84-a7ff-d7b9114ef1f2
# ╠═3b4822e4-5d27-4a1d-a093-2bb523608b49
# ╠═4c81a5a7-9e34-40b6-b9a5-e8138ad061bf
# ╠═662deff8-e425-4442-b524-7c656aeecbd6
# ╠═54c50afe-314b-46d4-a71d-ded8f028c1e7
# ╠═273b511b-330f-43dd-bf01-99bc980b9fd7
# ╠═f2e150ff-b33f-4697-bebd-790a6437fdf0
# ╠═448fdb1e-7a3c-4105-bd9b-f18d109986a4
