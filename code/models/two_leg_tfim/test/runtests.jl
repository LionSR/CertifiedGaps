# Fast unit tests for the two-leg TFIM ladder model (the stored certificates are checked by
# check.jl).
using Test, LinearAlgebra, CertifiedGaps, CertifiedGaps.RingModels
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum, float_ring_spectrum
const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "src", "Ladder.jl"))
using .Ladder
BLAS.set_num_threads(1)

@testset "Exact frame change and reflection positivity" begin
    W = rung_hadamard()
    @test W * W == I && W == transpose(W)
    for g in (2 // 1, 5 // 2, 3 // 1, 4 // 1)
        m = ladder(g)
        @test is_frame_change(physical_ladder(g), m, W)
        @test verify_rp(m, ladder_rp(m))
        # Physical frame: −h = Σ (Z⊗1) ⊗ (Z⊗1) + (1⊗Z) ⊗ (1⊗Z), as in the Lean model.
        Z, I2 = [1 0; 0 -1], [1 0; 0 1]
        @test -physical_terms(g).bond == kron(Z, I2, Z, I2) + kron(I2, Z, I2, Z)
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
        @test float_ring_spectrum(ladder(g), L).E ≈ E atol = 1e-11
        spec = verified_ring_spectrum(ladder(g), L)
        @test maximum(s.radius for s in spec) < 1e-9
        @test sum(length(s.centers) * s.multiplicity for s in spec) == N
    end
end

include(joinpath(ROOT, "src", "RitzSeeds.jl"))
using .RitzSeeds
using CertifiedGaps.RitzBounds: partial_ring_spectrum, log_partition_lower
using CertifiedGaps.VerifiedSpectra: log_partition_enclosure
using CertifiedGaps.GapCertificates: coupled_seed_bounds, seed_bounds

@testset "Coupled seeds from Ritz bounds (small rings)" begin
    # The doubled ring of 4 rungs with forced Lanczos in its larger blocks: the Ritz bounds
    # give lower bounds on log Z, and the coupled defects of the seed (2, β) are at least the
    # complete-spectrum ones, and close to them when every block keeps its low levels.
    g = 3 // 1
    m = ladder(g)
    full = verified_ring_spectrum(m, 4)
    part = partial_ring_spectrum(m, 4; howmany = ritz_levels, dense_max = 8)
    @test any(b.label["method"] == "ritz" for b in part)
    @test sum(b.label["dimension"] * b.multiplicity for b in part) == 4^4
    for β in (1 // 2, 2, 6)
        @test inf(log_partition_lower(part, β)) <= sup(log_partition_enclosure(full, β))
    end
    small = verified_ring_spectrum(m, 2)
    for β in (1, 2)
        coupled = coupled_seed_bounds(m, 2, β; small, large = part)
        complete = seed_bounds(m, 2, β; spectra = Dict(2 => small, 4 => full))
        @test coupled.eps_s >= complete.eps_s && coupled.eps_t >= complete.eps_t
        @test coupled.eps_s - complete.eps_s < 1e-3
    end
    @test ritz_levels(Dict("charge" => 0)) > ritz_levels(Dict("charge" => 1)) > 0
    @test all(p.n0 == RITZ_N0 && iseven(p.n0) for p in RITZ_POINTS)
    # The bridge of the 4-rung ring: an outward enclosure of its gap, far above the bounds.
    bridge = only(finite_bridge(g))
    lo, hi = (parse(Rational{BigInt}, bridge[k]) for k in ("gap_lower", "gap_upper"))
    @test bridge["L"] == 4 && 2.7 < lo <= hi < 2.73
end
