include(joinpath(@__DIR__, "../scripts/certify_trial.jl"))
using Test, Random

@testset "Exact spin-two intertwiners and multiplicity transfer" begin
    for spins in ([1], [1, 2], [1, 1, 2], [0, 1, 1, 2, 2, 3])
        t = ST.templates(spins)
        rng = MersenneTwister(61)
        x = randn(rng, length(t.pairs))
        A = ST.tensors(t, x)
        r = randn(rng, length(t.singlets))
        R = zeros(t.D, t.D)
        for (i, (a, b)) in enumerate(t.singlets), k in 1:t.dims[a]
            R[t.ranges[a][k], t.ranges[b][k]] = r[i]
        end
        direct = sum(a*R*transpose(a) for a in A)
        expected = zeros(t.D, t.D)
        wr = ST.singlet_transfer(t, x)*r
        for (i, (a, b)) in enumerate(t.singlets), k in 1:t.dims[a]
            expected[t.ranges[a][k], t.ranges[b][k]] = wr[i]
        end
        @test direct ≈ expected atol=1e-12
        @test ST.energy(t, x) ≈ ST.energy(t, 3x) atol=1e-11
        @test ST.energy(t, x) ≈ ST.energy(t, ST.canonicalize(t, x)) atol=1e-9
    end
    @test ST.energy(ST.templates([1]), [1.0]) ≈ -4.5 atol=1e-13
    @test_throws ArgumentError ST.intertwiner(0, 0)
end

@testset "Exact polynomial amplitudes and charge blocks" begin
    for spins in ([1], [1, 2])
        t = ST.templates(spins)
        P = ST.polynomial_tensors(t, big.(collect(1:length(t.pairs))))
        blocks = ST.exact_blocks(t, P)
        direct = ST.physical_integer_check(P)
        c = ST.exact_contract(blocks, 4)
        @test (c.V, 4c.U)==(direct.V, direct.E)
        # Full transfer powers independently check the charge assembly at another length.
        T = sum(ST.WEIGHTS[m]*kron(P[m], P[m]) for m in 1:5)
        @test sum(size(b.T, 1) for b in blocks)==t.D^2
        @test ST.exact_contract(blocks, 6).V == tr(T^6)
        @test_throws AssertionError ST.exact_blocks(t, [p*big(10)^30 for p in P])
    end
end

@testset "Saved finite spin-two energy certificates" begin
    source = joinpath(@__DIR__, "../results/trial-certificate.json")
    data = JSON3.read(read(source, String), Dict{String,Any})
    @test data["status"]=="finite-energy-certificates-passed"
    @test data["uniform_gap_certified"]==false
    for (name, digest) in data["source_sha256"]
        @test bytes2hex(sha256(read(joinpath(@__DIR__, "../src", name))))==digest
    end
    t = ST.templates(Int.(data["virtual_spins"]))
    parameters = parse.(BigInt, data["integer_parameters"])
    P = ST.polynomial_tensors(t, parameters)
    direct = ST.physical_integer_check(P)
    @test direct.V>0
    control = itensors_check(P, direct)
    @test control["absolute_energy_difference"]<1e-10
    @test control["total_raising_residual"]<1e-6
    @test [c["length"] for c in data["certificates"]]==[128, 256, 512]
    for c in data["certificates"]
        U, V = parse(BigInt, c["U"]), parse(BigInt, c["V"])
        @test V>0 && 10000U+47595V<0
        @test -47595V-10000U == parse(BigInt, c["positive_margin_minus47595V_minus10000U"])
        @test fld(big(10)^12*U, V) == parse(BigInt, c["energy_per_bond_floor_1e12"])
        @test sum(parse(BigInt, b["U"]) for b in c["blocks"])==U
        @test sum(parse(BigInt, b["V"]) for b in c["blocks"])==V
    end
    x = Float64.(parameters) .* t.norms
    @test abs(ST.energy(t, x)-last(data["certificates"])["energy_per_bond_display"])<1e-9
end
