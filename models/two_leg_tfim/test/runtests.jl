# Fast tests and certificate checks for the two-leg TFIM ladder (run by check.jl).
using Test, LinearAlgebra, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.RingModels
using CertifiedGaps.GapCertificates: check_certificate
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "src", "Ladder.jl"))
using .Ladder
BLAS.set_num_threads(1)

certificate_paths() = filter(p -> occursin(r"^certificate_g.*\.json$", basename(p)),
    readdir(joinpath(ROOT, "results"); join = true))

# The seed defects hardcoded as the first entries of `s` and `t` in
# lean/CertifiedGaps/Models/TwoLegIsing.lean, per field.
const LEAN_SEEDS = Dict(3 // 1 => (s = 1125 // 10^4, t = 5162 // 10^5))

@testset "Exact frame change and reflection positivity" begin
    W = rung_hadamard()
    @test W * W == I && W == transpose(W)
    for g in (2 // 1, 5 // 2, 3 // 1, 4 // 1)
        p, r = physical_terms(g), rotated_terms(g)
        @test W * p.onsite * W == r.onsite
        @test kron(W, W) * p.bond * kron(W, W) == r.bond
        m = ladder(g)
        @test verify_rp(m, ladder_rp(m))
        # Physical frame: −h = Σ (Z⊗1) ⊗ (Z⊗1) + (1⊗Z) ⊗ (1⊗Z), as in the Lean model.
        Z, I2 = [1 0; 0 -1], [1 0; 0 1]
        @test -p.bond == kron(Z, I2, Z, I2) + kron(I2, Z, I2, Z)
    end
end

@testset "Core blocks against an independent full Hamiltonian" begin
    # Rotated frame, rung state index 2a + b + 1 with a, b ∈ {0 = up, 1 = down}; spin (j, leg) is
    # bit 2j + leg of the ring state s, leg a the more significant bit of each rung.
    g = 3
    for L in 2:4
        N = 4^L; H = zeros(N, N)
        bit(j, leg) = 2j + (1 - leg)
        for s in 0:(N-1)
            H[s+1, s+1] = -g * (2L - 2count_ones(s))
            for j in 0:(L-1)
                H[(s ⊻ (1 << bit(j, 0)) ⊻ (1 << bit(j, 1)))+1, s+1] -= 1
                for leg in 0:1
                    H[(s ⊻ (1 << bit(j, leg)) ⊻ (1 << bit(mod(j + 1, L), leg)))+1, s+1] -= 1
                end
            end
        end
        E = eigvals(Symmetric(H))
        @test ladder_energies(g, L) ≈ E atol = 1e-11
        spec = verified_ring_spectrum(ladder(g), L)
        @test maximum(s.radius for s in spec) < 1e-9
        @test sum(length(s.centers) * s.multiplicity for s in spec) == N
    end
end

@testset "Stored certificates" begin
    paths = certificate_paths()
    @test !isempty(paths)
    for path in paths
        c = read_exact_json(path)
        g = parse(Rational{BigInt}, c["point"]["g"])
        point = only(filter(p -> p.g == g, POINTS))
        # The stored model and RP data are those of `ladder(g)`, with the listed seed.
        stored, model = decode_model(c["model"]), ladder(g)
        @test (stored.d, stored.bond, stored.onsite, stored.charge, stored.modulus,
            stored.symmetries) ==
              (model.d, model.bond, model.onsite, model.charge, model.modulus,
            model.symmetries)
        @test parse(Int, c["seed"]["n0"]) == point.n0
        @test parse(Rational{BigInt}, c["seed"]["beta0"]) == point.β0
        # Exact RP identity and the exact chain from the stored seeds (throws on mismatch).
        @test check_certificate(c; spectra = false)["gap_lower_bound"] == c["gap_lower_bound"]
        εs = parse(Rational{BigInt}, c["seed"]["eps_s"])
        εt = parse(Rational{BigInt}, c["seed"]["eps_t"])
        @test 0 <= εs < 1 && 0 <= εt < 1
        @test certified_gap_bound(εs, εt, point.β0).bound == c["interval_bound_decimal"]
        # The Lean theorem's seed hypotheses are implied by the certified defects.
        @test εs <= LEAN_SEEDS[g].s && εt <= LEAN_SEEDS[g].t
        for (L, s) in c["spectra"]
            @test parse(Int, s["dimension"]) == 4^parse(Int, L)
        end
        @test verify_model_metadata(c["metadata"], ROOT; PROVENANCE...) == 3
    end
end
