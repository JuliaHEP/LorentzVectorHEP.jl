using LorentzVectorHEP
using Test

@testset "calculations" begin
    v1 = LorentzVectorCyl(1761.65,-2.30322,-2.5127,0.105652)
    v2 = LorentzVectorCyl(115.906,-2.28564,-2.50781,0.105713)
    @test energy(v1) ≈ 8901.870789524375 atol=1e-6
    @test px(v1) ≈ -1424.610065192358 atol=1e-6
    @test py(v1) ≈ -1036.2899616674022 atol=1e-6
    @test pz(v1) ≈ -8725.817601790963 atol=1e-6

    @test isapprox((v1+v2).mass, 8.25741602000877, atol=1e-6)
    @test isapprox(fast_mass(v1,v2), 8.25741602000877, atol=1e-6)

    v1 = LorentzVectorCyl(43.71242f0, 1.4733887f0, 1.6855469f0, 0.10571289f0)
    v2 = LorentzVectorCyl(36.994347f0, 0.38684082f0, -1.3935547f0, 0.10571289f0)
    @test (v1+v2).mass == 92.55651f0
    @test fast_mass(v1,v2) == 92.55651f0
    @test isapprox(deltar(v1,v2), 3.265188f0, atol=1e-6)
    @test isapprox(deltaphi(v1,v2), -3.0791016f0, atol=1e-6)

    v3 = v1*5
    @test v3.pt == 5*v1.pt
    @test v3.mass == 5*v1.mass
    @test v3.eta == v1.eta
    @test v3.phi == v1.phi
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
    v2 = LorentzVectorCyl(LorentzVector(v1))
    @test v1.pt ≈ v2.pt
    @test v1.eta ≈ v2.eta
    @test v1.phi ≈ v2.phi
    @test v1.mass ≈ v2.mass atol=1e-6
    for func in (px, py, pz, energy, pt, eta, phi, mass)
        func(v1) ≈ func(LorentzVector(v1))
    end
end
