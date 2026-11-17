using Test, LinearAlgebra, IntervalArithmetic
using CertifiedGaps
using CertifiedGaps.ExactArithmetic: frac, parse_rational, round_up
using CertifiedGaps.Results: read_result, write_result
using CertifiedGaps.VerifiedSpectra
include(joinpath(@__DIR__, "..", "src", "Certificate.jl"))
using .Certificate
using .Certificate.SquareIsing, .Certificate.CubeRecurrence

const ROOT = dirname(@__DIR__)
BLAS.set_num_threads(1)

# Independently enumerate spins and nearest-neighbor edges, without row operators.
function direct_hamiltonian(n, m, g)
    N = n*m
    H = zeros(Int, 2^N, 2^N)
    bit(x, y) = mod(x, n)*m + mod(y, m)
    for s in 0:2^N-1
        H[s+1, s+1] = -g*(N-2count_ones(s))
        for x in 0:n-1, y in 0:m-1, (dx, dy) in ((1, 0), (0, 1))
            t = xor(s, (1 << bit(x, y)), (1 << bit(x+dx, y+dy)))
            H[t+1, s+1] -= 1
        end
    end
    return H
end

@testset "independent torus geometry and complete small spectra" begin
    for (n, m, sum_weights_squared) in ((2, 2, 16), (4, 2, 24))
        H = direct_hamiltonian(n, m, 7)
        @test H == transpose(H)
        @test tr(H) == 0
        @test sum(abs2, H) == 2^(n*m)*(n*m*7^2 + sum_weights_squared)
        E = eigvals(Symmetric(Float64.(H)))
        @test float_ring_spectrum(row_model(m, 7), n).E ≈ E atol = 1e-10
        spectrum = verified_ring_spectrum(row_model(m, 7), n)
        @test sum(length(s.centers)*s.multiplicity for s in spectrum) == 2^(n*m)
        @test maximum(s.radius for s in spectrum) < 1e-8
    end
    # Interchange the spatial directions; the two repeated-bond conventions must agree.
    @test eigvals(Symmetric(Float64.(direct_hamiltonian(2, 4, 7)))) ≈
        float_ring_spectrum(row_model(2, 7), 4).E atol = 1e-10
    @test row_rp(4, 7).weights == ones(Int, 4)
end

@testset "cube arithmetic and finite windows" begin
    small = (frac(1, 1000), frac(1, 1000), frac(1, 1000))
    @test maximum(cube_step(small)) <= 25maximum(small)^2
    chain = cube_chain(small, frac(1))
    @test frac(3688, 1000) < parse_rational(chain["gap_lower_bound"]) < frac(3689, 1000)
    asymmetric = (frac(1, 100), frac(2, 100), frac(1, 10000))
    # Guard the staggered order independently of the shared cube-product helper.
    for (x, y, t) in (small, asymmetric), grid in (big(10)^40, big(10)^80)
        update(z) = z^2 / (2(1-z)^2)
        a = round_up(update(1 - (1-x)^4 * (1-y)^2 * (1-t)), grid)
        b = round_up(update(1 - (1-y)^4 * (1-t)^2 * (1-a)), grid)
        c = round_up(update(1 - (1-t)^4 * (1-a)^2 * (1-b)), grid)
        @test cube_step((x, y, t); grid) == (a, b, c)
    end
    chain = cube_chain(asymmetric, frac(1, 5))
    @test length(chain["rows"]) > 1
    @test parse_rational(chain["gap_lower_bound"]) > 0
    @test_throws ErrorException cube_chain((frac(1, 2), frac(1, 2), frac(1, 2)), frac(1))
    @test_throws ErrorException cube_chain(small, frac(0))
    @test_throws ErrorException cube_step((frac(-1), frac(0), frac(0)))
    refined = cube_chain(asymmetric, frac(1, 5); min_steps = 8, grid = big(10)^400)
    @test parse_rational(refined["gap_lower_bound"]) > parse_rational(chain["gap_lower_bound"])
    @test length(refined["rows"]) == 9
    @test_throws ErrorException cube_chain(small, frac(1); min_steps = 13)
    @test_throws ErrorException cube_step(small; grid = 0)
end

@testset "direct-temperature prefix and parameter domains" begin
    ex, et = frac(1, 10), frac(3, 100)
    grid = big(10)^12
    f(z) = z^2 / (2*(1-z)^2)
    a = round_up(f(ex), grid)
    b = round_up(f(1-(1-ex)^2*(1-a)), grid)
    c = round_up(f(1-(1-et)^4*(1-a)^2*(1-b)), grid)
    @test direct_temperature_prefix(ex, et; grid) == (a, b, c)
    @test_throws ErrorException direct_temperature_prefix(-ex, et)
    @test_throws ErrorException direct_temperature_prefix(ex, 1)
    @test_throws ErrorException direct_temperature_prefix(frac(9, 10), et)
    @test_throws ErrorException direct_temperature_prefix(ex, et; grid = 0)
    @test_throws ErrorException cube_chain((ex, et, et), 1; min_even_side = 3)
    @test_throws ErrorException Certificate.direct_seed_chain([ex, et], 0)
    @test_throws ErrorException Certificate.parameter_interval_record(15, 1, [ex, et]; radius = 0)
    @test_throws ErrorException Certificate.parameter_interval_record(15, 1, [ex, et]; radius = 15)
    @test_throws ErrorException Certificate.parameter_interval_record(15, 1, [ex, et]; radius = 1)
end

# Check the defining inequalities without calling cube_step, cube_chain, or the prefix helper.
function audit_chain(chain, beta, side)
    rows = chain["rows"]
    @test chain["min_even_side"] == string(side)
    @test length(rows) == 13
    f(z) = z^2 / (2*(1-z)^2)
    claimed = parse_rational(chain["gap_lower_bound"])
    for (i, row) in enumerate(rows)
        x, y, t = parse_rational.(row["defects"])
        @test all(z -> 0 < z < 1, (x, y, t))
        @test parse_rational(row["beta"]) == beta*big(2)^(i-1)
        if i == length(rows)
            r = max(x, y, t)
            @test r <= frac(1, 1000)
            lower, _ = CertifiedGaps.ExactArithmetic.log_enclosure(1/(25r))
            @test claimed <= lower/(beta*big(2)^(i-1))
        else
            Q = (1-x)^4*(1-y)^2*(1-t)
            @test Q > frac(1, 2)
            @test parse_rational(row["cube_lower"]) == Q
            lower, _ = CertifiedGaps.ExactArithmetic.log_enclosure(Q/(1-Q))
            @test claimed <= lower/(beta*big(2)^(i-1))
            a, b, c = parse_rational.(rows[i+1]["defects"])
            @test a >= f(1-Q)
            @test b >= f(1-(1-y)^4*(1-t)^2*(1-a))
            @test c >= f(1-(1-t)^4*(1-a)^2*(1-b))
        end
    end
end

@testset "stored certificates, independent inequalities, and spectral moments" begin
    minimum_gaps = Dict(15 => frac(74208, 100000), 16 => frac(17235, 10000),
        18 => frac(33912, 10000), 20 => frac(49771, 10000), 24 => frac(81419, 10000))
    for p in POINTS
        record = Certificate.read_certificate(ROOT, p.g)
        @test check(ROOT, p.g) == record
        @test uniform_bound(record) >= minimum_gaps[p.g]
        @test parse_rational(record["finite_volume"]["4x4"]["gap"][1]) > uniform_bound(record)
        direct = record["direct_temperature"]
        ex, et = parse_rational.(direct["input"]["defects"])
        a, b, c = parse_rational.(direct["seed"]["defects"])
        f(z) = z^2 / (2*(1-z)^2)
        @test a >= f(ex)
        @test b >= f(1-(1-ex)^2*(1-a))
        @test c >= f(1-(1-et)^4*(1-a)^2*(1-b))
        @test direct["seed"]["n0"] == direct["seed"]["m0"] == "4"
        audit_chain(direct["chain"], 2p.beta, 4)
        region = record["parameter_interval"]
        lo, hi = parse_rational.(region["field_interval"])
        r = (p.g-lo)/(p.g*lo)
        loss = 64*(2p.beta*p.g)*r
        @test p.g-lo == hi-p.g == frac(1, 1000)
        @test parse_rational(region["log_purity_loss"]) == loss
        @test parse_rational.(region["input_defects"]) == [1-(1-x)*(1-loss) for x in (ex, et)]
        @test parse_rational(region["gap_lower_bound"]) ==
            lo/p.g * parse_rational(region["chain"]["gap_lower_bound"])
        audit_chain(region["chain"], 2p.beta, 4)
        if p.g in (18, 20, 24)
            @test record["legacy"]["chain"]["min_even_side"] == "2"
            @test parse_rational(record["legacy"]["refinement"]["chain"]["gap_lower_bound"]) <
                uniform_bound(record)
        else
            @test isnothing(record["legacy"])
        end
        for (n, m, sw2) in ((4, 2, 24), (4, 4, 32))
            file = record["spectral_record"]["spectra"]["$(n)x$(m)"]["file"]
            spectrum = read_spectrum(joinpath(ROOT, "results", file))
            trace1, trace2 = interval(0.0), interval(0.0)
            for block in spectrum, center in block.centers
                radius = sup(interval(block.radius))
                eigenvalue = interval(center) + interval(-radius, radius)
                trace1 += block.multiplicity * eigenvalue
                trace2 += block.multiplicity * eigenvalue^2
            end
            @test inf(trace1) <= 0 <= sup(trace1)
            @test inf(trace2) <= 2^(n*m)*(n*m*p.g^2+sw2) <= sup(trace2)
        end
    end
end

@testset "reject altered scientific claims and evidence" begin
    original = Certificate.read_certificate(ROOT, 15)
    mutations = [
        r -> r["g"] = "16",
        r -> r["unexpected"] = "claim",
        r -> r["metadata"]["source_sha256"]["models/square_tfim/src/Certificate.jl"] = "bad",
        r -> r["spectral_record"]["models"]["4"] = Dict("model" => "wrong"),
        r -> r["direct_temperature"]["chain"]["min_even_side"] = "2",
        r -> r["direct_temperature"]["chain"]["rows"][1]["cube_lower"] = "1",
        r -> r["scope"] = "unique ground state; every even square torus of side L >= 2",
        r -> r["direct_temperature"]["input"]["beta"] = 0.235,
        r -> r["direct_temperature"]["seed"]["defects"][1] = "1//100000000000",
        r -> r["direct_temperature"]["chain"]["gap_lower_bound"] = "1000000",
        r -> r["parameter_interval"]["field_interval"][1] = "14",
        r -> r["parameter_interval"]["log_purity_loss"] = "0",
        r -> r["parameter_interval"]["gap_lower_bound"] = "1000000",
        r -> r["finite_volume"]["4x4"]["ground_energy"][1] = "0",
        r -> r["spectral_record"]["g"] = "16",
        r -> r["spectral_record"]["spectra"]["4x4"]["sha256"] = "bad",
        r -> r["spectral_record"]["model_source_sha256"] = "bad",
    ]
    for mutate in mutations
        changed = deepcopy(original)
        mutate(changed)
        @test_throws ErrorException Certificate.check_record(ROOT, 15, changed)
    end
    changed = deepcopy(Certificate.read_certificate(ROOT, 18))
    changed["legacy"]["refinement"]["chain"]["gap_lower_bound"] = "1000000"
    @test_throws ErrorException Certificate.check_record(ROOT, 18, changed)
    # The exact reader must not send the long, already-quoted rationals through Float64 or regex.
    mktempdir() do temporary
        mkpath(joinpath(temporary, "results"))
        write_result(certificate_path(temporary, 15), original)
        @test Certificate.read_certificate(temporary, 15) == original
    end
end
