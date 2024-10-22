using LorentzVectorBase
using LorentzVectorHEP
using Test

include("patch_test.jl")

@testset "calculations" begin
    v1 = LorentzVectorCyl(1761.65,-2.30322,-2.5127,0.105652)
    v2 = LorentzVectorCyl(115.906,-2.28564,-2.50781,0.105713)
    @test LorentzVectorBase.energy(v1) ≈ 8901.870789524375 atol=1e-6
    @test LorentzVectorBase.px(v1) ≈ -1424.610065192358 atol=1e-6
    @test LorentzVectorBase.py(v1) ≈ -1036.2899616674022 atol=1e-6
    @test LorentzVectorBase.pz(v1) ≈ -8725.817601790963 atol=1e-6
    @test LorentzVectorBase.rapidity(v1) ≈ -2.3032199982371715 atol=1e-6

    @test LorentzVectorBase.pt2(v1) ≈ 3.1034107225e6 atol=1e-6
    @test LorentzVectorBase.mass2(v1) ≈ 0.011162345103999998 atol=1e-6

    @test isapprox((v1+v2).mass, 8.25741602000877, atol=1e-6)
    @test isapprox(fast_mass(v1,v2), 8.25741602000877, atol=1e-6)

    v1 = LorentzVectorCyl(43.71242f0, 1.4733887f0, 1.6855469f0, 0.10571289f0)
    v2 = LorentzVectorCyl(36.994347f0, 0.38684082f0, -1.3935547f0, 0.10571289f0)
    @test (v1+v2).mass == 92.55651f0
    @test fast_mass(v1,v2) ≈ 92.55651f0
    @test isapprox(deltar(v1,v2), 3.265188f0, atol=1e-6)
    @test isapprox(deltaphi(v1,v2), -3.0791016f0, atol=1e-6)

    v3 = v1*5
    @test v3.pt == 5*v1.pt
    @test v3.mass == 5*v1.mass
    @test v3.eta == v1.eta
    @test v3.phi == v1.phi

    vcart1 = old_LorentzVector(10.0, -2.3, 4.5, 0.23)
    @test LorentzVectorBase.rapidity(vcart1) ≈ 0.02300405695442185 atol=1e-9
    @test LorentzVectorBase.eta(vcart1) ≈ 0.045495409709778126 atol=1e-9
    @test LorentzVectorBase.phi(vcart1) ≈ 2.0432932623119604 atol=1e-9

    vcart2 = old_LorentzVector(10.0, 2.7, -4.1, -0.21)
    @test LorentzVectorBase.rapidity(vcart2) ≈ -0.021003087817077763 atol=1e-9
    @test LorentzVectorBase.eta(vcart2) ≈ -0.04276400891568771 atol=1e-9
    @test LorentzVectorBase.phi(vcart2) ≈ -0.9884433806509134 atol=1e-9
    @test LorentzVectorHEP.phi02pi(vcart2) ≈ 5.294741926528673 atol=1e-9

    @test vcart1 + vcart2 ≈ old_LorentzVector(vcart1.t + vcart2.t, vcart1.x + vcart2.x, vcart1.y + vcart2.y, vcart1.z + vcart2.z)
    @test vcart1 + vcart2 == vcart2 + vcart1
    @test -vcart1 == old_LorentzVector(-(vcart1.t), -(vcart1.x), -(vcart1.y), -(vcart1.z))
    @test vcart1 - vcart2 ≈ old_LorentzVector(vcart1.t - vcart2.t, vcart1.x - vcart2.x, vcart1.y - vcart2.y, vcart1.z - vcart2.z)
    @test vcart1 - vcart2 == -(vcart2 - vcart1)
    
    c = rand()
    @test c*vcart1 ≈ old_LorentzVector(c*vcart1.t, c*vcart1.x, c*vcart1.y, c*vcart1.z)
    @test -c*vcart1 ≈ old_LorentzVector(-c*vcart1.t, -c*vcart1.x, -c*vcart1.y, -c*vcart1.z)
    @test vcart1*c == c*vcart1

    if c == 0
        c += 0.1
    end

    @test vcart1 / c ≈ vcart1 * (1/c)
    

    @test deltaeta(vcart1, vcart2) ≈ 0.08825941862546584 atol=1e-9
    @test deltaphi(vcart1, vcart2) ≈ 3.0317366429628736 atol=1e-9

    vcart3 = old_LorentzVector(66.0, 0.0, 0.0, 66.0)
    @test_broken rapidity(vcart3) ≈ 100066.0 atol=1e-9

    vcart4 = old_LorentzVector(4.4, 8.1, 2.2, 3.3)
    @test LorentzVectorBase.mass(vcart4) ≈ -7.872737770305829 atol=1e-9
end

@testset "summing" begin
    v1 = LorentzVectorCyl(1761.65,-2.30322,-2.5127,0.105652)
    v2 = LorentzVectorCyl(115.906,-2.28564,-2.50781,0.105713)
    v3 = LorentzVectorCyl(43.71242f0, 1.4733887f0, 1.6855469f0, 0.10571289f0)
    v4 = LorentzVectorCyl(36.994347f0, 0.38684082f0, -1.3935547f0, 0.10571289f0)
    vs = [v1, v2, v3, v4]
    @test sum(vs).mass ≈ 2153.511000993
    @test sum(LorentzVectorCyl[]).mass ≈ 0
end

@testset "broadcasting" begin
    pts = [1761.65,115.906,43.712420,36.994347]
    etas = [-2.30322,-2.28564,1.4733887,0.38684082]
    phis = [-2.5127,-2.50781,1.6855469,-1.3935547]
    mass = 0.105652
    vs = LorentzVectorCyl.(pts, etas, phis, mass)
    @test all([v.mass for v in vs] .== mass)
    @test fast_mass.(vs[1], vs[2:end]) == fast_mass.(Ref(vs[1]), vs[2:end])
end

@testset "conversions" begin
    v1 = LorentzVectorCyl(1761.65,-2.30322,-2.5127,0.105652)
    v2 = LorentzVectorCyl(old_LorentzVector(v1))
    @test v1.pt ≈ v2.pt
    @test v1.eta ≈ v2.eta
    @test v1.phi ≈ v2.phi
    @test v1.mass ≈ v2.mass atol=1e-6
    for func in (LorentzVectorBase.px, LorentzVectorBase.py, LorentzVectorBase.pz, LorentzVectorBase.energy, LorentzVectorBase.pt, LorentzVectorBase.eta, LorentzVectorBase.phi, LorentzVectorBase.mass)
        func(v1) ≈ func(old_LorentzVector(v1))
    end
end

@testset "fromPtEtaPhiE" begin
    v1 = LorentzVectorCyl(1761.65,-2.30322,-2.5127,0.105652) 
    v2 = fromPtEtaPhiE(v1.pt, v1.eta, v1.phi, LorentzVectorBase.energy(v1))
    @test v1.mass ≈ v2.mass atol=1e-6
 end
