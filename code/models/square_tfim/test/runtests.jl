using Test, LinearAlgebra, IntervalArithmetic
using CertifiedGaps.ExactArithmetic: frac, parse_rational
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

@testset "stored certificate and rejection of altered claims" begin
    for path in sort(filter(p -> occursin(r"certificate_g[0-9]+\.json\.gz$", p),
            readdir(joinpath(ROOT, "results"); join = true)))
        record = read_result(path)
        g = parse(Int, record["g"])
        @test check(ROOT, g)["chain"] == record["chain"]
        @test uniform_bound(record) > parse_rational(record["chain"]["gap_lower_bound"])
        @test parse_rational(record["finite_volume"]["4x4"]["gap"][1]) > uniform_bound(record)
        name = record["spectral_record"]["spectra"]["4x4"]["file"]
        spectrum = read_spectrum(joinpath(ROOT, "results", name))
        trace1, trace2 = interval(0.0), interval(0.0)
        for block in spectrum, center in block.centers
            r = sup(interval(block.radius))
            eigenvalue = interval(center) + interval(-r, r)
            trace1 += block.multiplicity * eigenvalue
            trace2 += block.multiplicity * eigenvalue^2
        end
        @test inf(trace1) <= 0 <= sup(trace1)
        @test inf(trace2) <= 2^16*(16g^2+32) <= sup(trace2)
        mktempdir() do temporary
            modelroot = joinpath(temporary, "code", "models", "square_tfim")
            mkpath(dirname(modelroot))
            cp(ROOT, modelroot; follow_symlinks = true)
            # Provenance compares with an actual core environment, not a fabricated hash set.
            symlink(normpath(joinpath(ROOT, "..", "..", "julia")),
                joinpath(temporary, "code", "julia"))
            @test check(modelroot, g)["chain"] == record["chain"]
            altered = deepcopy(record)
            altered["chain"]["gap_lower_bound"] = "1000000"
            write_result(certificate_path(modelroot, g), altered)
            @test_throws ErrorException check(modelroot, g)
            altered = deepcopy(record)
            altered["refinement"]["chain"]["gap_lower_bound"] = "1000000"
            write_result(certificate_path(modelroot, g), altered)
            @test_throws ErrorException check(modelroot, g)
            altered = deepcopy(record)
            altered["finite_volume"]["4x4"]["ground_energy"][1] = "0"
            write_result(certificate_path(modelroot, g), altered)
            @test_throws ErrorException check(modelroot, g)
            altered = deepcopy(record)
            altered["seed"]["defects"][1] = "1//100000000000"
            write_result(certificate_path(modelroot, g), altered)
            @test_throws ErrorException check(modelroot, g)
            altered = deepcopy(record)
            altered["spectral_record"]["spectra"]["4x4"]["sha256"] = "bad"
            write_result(certificate_path(modelroot, g), altered)
            @test_throws ErrorException check(modelroot, g)
        end
    end
end
