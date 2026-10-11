include(joinpath(@__DIR__, "../scripts/one_trace_target.jl"))
using Test

@testset "Conditional one-trace lifting" begin
    for (num, h, r) in ((102, 2, 2), (104, 2, 2), (106, 2, 3), (108, 2, 3))
        U = OT.frac(num, 100)
        target = OT.target(U, 128, OT.frac(32); h, r)
        @test !isnothing(target)
        @test !target["spin2_thermal_input_certified"]
        @test target["trial_length_required"] == 512
        @test target["conditional_gap_proof"]["gap_strict_lower_decimal"] > 0
        lower_s = OT.parse_rational(target["spatial_purity_lower"])
        lower_t = OT.parse_rational(target["temporal_purity_lower"])
        D, Q = U^r, U^(2h)
        for count in (1, 3, 11)
            # Positive moment vectors satisfying total <= D and 2h-moment > 1.
            weights = vcat(OT.frac(1), fill((D-1)/count, count))
            @test sum(weights) == D && sum(w^(2h) for w in weights) > 1
            spatial = sum(w^(2h) for w in weights)/sum(w^h for w in weights)^2
            @test spatial >= lower_s
            # Gibbs ratios at beta whose sum saturates the allowed excited mass.
            excitations = fill((Q-1)/count, count)
            temporal = (1+sum(e^(2r) for e in excitations))/(1+sum(e^r for e in excitations))^2
            @test temporal >= lower_t
        end
    end
    @test isnothing(OT.target(OT.frac(3, 2), 128, OT.frac(32); h = 2, r = 2))
    trial = JSON3.read(
        read(joinpath(@__DIR__, "../results/trial-137-certificate.json"), String),
        Dict{String,Any},
    )
    finite = only(filter(c->c["length"]==128, trial["certificates"]))["energy_per_bond_interval"]
    upper = parse(BigInt, finite["upper_numerator"])//parse(BigInt, finite["denominator"])
    @test 1-32*128*(OT.frac(238061, 50000)+upper) > OT.frac(11967, 10000)
    current = JSON3.read(
        read(joinpath(@__DIR__, "../results/trial-295-symmetric-certificate.json"), String),
        Dict{String,Any},
    )
    @test current["complete"] && !current["uniform_gap_certified"]
    shift = OT.frac(47612463, 10000000)
    for row in current["certificates"]
        interval = row["energy_per_bond_interval"]
        bound = parse(BigInt, interval["upper_numerator"])//parse(BigInt, interval["denominator"])
        @test bound < -shift
        if row["length"] == 256
            @test 1-64*256*(shift+bound) < OT.frac(102, 100)
        end
    end
    candidate = OT.target(OT.frac(108, 100), 256, OT.frac(64); h = 2, r = 3)
    @test candidate["trial_length_required"] == 1024
    @test candidate["conditional_gap_proof"]["gap_strict_lower_decimal"] > 0.003019
    @test !candidate["spin2_thermal_input_certified"]
end
