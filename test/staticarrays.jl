using LorentzVectorHEP, StaticArrays
using Test

@testset "static array conversions" begin
    for LV ∈ (LorentzVector, LorentzVectorCyl)
        v1 = LV(5.0, 0.1, 0.2, 0.3)
    
        v2 = SArray{Tuple{4}}(v1)
        @test v2 == [5.0, 0.1, 0.2, 0.3]
        v3 = StaticArray{Tuple{4}}(v1)
        @test v3 == [5.0, 0.1, 0.2, 0.3]
        v4 = SArray(v1)
        @test v4 == [5.0, 0.1, 0.2, 0.3]
        v5 = SVector{4}(v1)
        @test v5 == [5.0, 0.1, 0.2, 0.3]
    end

    v = LorentzVector(SVector{4}(5.0,0.1,0.2,0.3))
    @test energy(v) == 5.0
    @test px(v) == 0.1
    @test py(v) == 0.2
    @test pz(v) == 0.3

    v = LorentzVectorCyl(SVector{4}(5.0,0.1,0.2,0.3))
    @test pt(v) == 5.0
    @test eta(v) == 0.1
    @test phi(v) == 0.2
    @test mass(v) == 0.3
end
