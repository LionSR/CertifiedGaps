# Fast checks of the five core certificates: the exact sublattice rotation to the RP frame, the
# stored model, RP data and blocking against their definitions, the exact rational chain
# (`check_certificate(cert; spectra = false)`), agreement with the seeds hardcoded in the Lean
# theorem, model-local source hashes, and the compact summary. `check.jl --full` also re-verifies
# every seed spectrum.
using Test
include(joinpath(@__DIR__, "..", "src", "common.jl"))
using CertifiedGaps.ExactDiagonalization: spectrum
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum
const REPO = normpath(joinpath(ROOT, "..", ".."))

# The two-site operator B (on kron(site i, site j)) on a ring of L sites of dimension d, exactly;
# site 1 is the most significant tensor factor.
function two_site(B, i, j, L, d)
    H = zeros(eltype(B), d^L, d^L)
    for c in 0:(d^L-1)
        state = reverse(digits(c; base = d, pad = L))
        col = state[i] * d + state[j] + 1
        for row in 1:d^2
            iszero(B[row, col]) && continue
            new = copy(state)
            new[i], new[j] = divrem(row - 1, d)
            H[foldl((a, x) -> a * d + x, new; init = 0)+1, c+1] += B[row, col]
        end
    end
    return H
end
# Σ_j B_{j,j+1} + Σ_j o(j)_j on the ring of L sites (for L = 2 both bonds (1,2) and (2,1)).
ring_matrix(B, o, L) = sum(two_site(B, j, mod1(j + 1, L), L, 2) +
                           two_site(kron(o(j), I2), j, mod1(j + 1, L), L, 2) for j in 1:L)

@testset "exact sublattice rotation to the RP frame" begin
    @test U * U' == I2 && U * X * U' == -X && U * J * U' == J && U * Z * U' == -Z
    @test kron(U, I2) * physical_bond() * kron(U, I2)' == rotated_bond()
    @test kron(I2, U) * physical_bond() * kron(I2, U)' == rotated_bond()
    for h in FIELDS, L in (2, 4, 6)
        V = reduce(kron, [isodd(j) ? U : I2 for j in 1:L])
        # Physical H = Σ S·S − h Σ_j (−1)^j Sᶻ_j on sites j = 1, …, L.
        physical = ring_matrix(physical_bond(), j -> -h * (-1)^j * Z, L)
        m = site_model(h)
        rotated = ring_matrix(real(m.bond), j -> real(m.onsite), L)
        @test V * physical * V' == rotated
    end
end

@testset "RP data and the two-site cell charge" begin
    for h in FIELDS
        m = site_model(h)
        @test verify_rp(m, site_rp(h)) && isreal_model(m)
        # In cells (a, b) the charge Sᶻ_a − Sᶻ_b is minus the physical Sᶻ (U flips site a).
        @test BLOCKING.charge == [Int(Z[a, a] - Z[b, b]) for a in 1:2 for b in 1:2]
    end
end

@testset "verified spectra against physical-frame ED" begin
    for L in (4, 6)
        spec = verified_ring_spectrum(site_model(3 // 4), L)
        levels = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spec]))
        ref = spectrum(CertifiedGaps.Models.staggered_heisenberg(3 // 4).chain, L;
            dense_max = 10^6).E
        @test maximum(abs, levels - ref) < 1e-12
    end
end

# Seeds hardcoded in lean/CertifiedGaps/Models/StaggeredHeisenberg.lean: per point ⟨h, β₀, …⟩
# the first entries s[0], t[0] of the chain lists, which the Lean theorem takes as seed defects.
function lean_points()
    text = read(joinpath(REPO, "lean", "CertifiedGaps", "Models", "StaggeredHeisenberg.lean"),
        String)
    q(s) = (p = parse.(BigInt, split(s, "/")); length(p) == 1 ? p[1] // 1 : p[1] // p[2])
    sci(m, e) = parse(BigInt, m) // big(10)^parse(Int, e)
    pattern = r"⟨([0-9/ ]+), ([0-9/ ]+), \{\s*s := \[(\d+)e-(\d+),.*?t := \[(\d+)e-(\d+),"s
    return [(h = q(replace(m[1], " " => "")), β = q(replace(m[2], " " => "")),
        s = sci(m[3], m[4]), t = sci(m[5], m[6])) for m in eachmatch(pattern, text)]
end

@testset "certificates (exact parts)" begin
    lean = lean_points()
    @test length(lean) == length(SEEDS)
    for (h, β0) in SEEDS
        cert = read_exact_json(certificate_path(h))
        @test parse_rational(cert["staggered_field"]) == h
        # The stored model and RP data are the intended ones, not merely self-consistent.
        model = decode_model(cert["model"])
        reference = site_model(h)
        @test model.d == 2 && model.bond == reference.bond && model.onsite == reference.onsite
        rp = decode_rp(cert["reflection_positivity"])
        @test rp.factors == site_rp(h).factors && rp.reflection === nothing
        @test parse(Int, cert["seed"]["n0"]) == N0 && parse_rational(cert["seed"]["beta0"]) == β0
        b = cert["blocking"]
        cell = decode_model(b["model"])
        @test parse(Int, b["cell"]) == BLOCKING.cell && cell.d == 4
        @test cell.charge == BLOCKING.charge && cell.modulus == BLOCKING.modulus
        recomputed = check_certificate(cert; spectra = false)
        @test recomputed["gap_lower_bound"] == cert["gap_lower_bound"]
        @test parse_rational(cert["gap_lower_bound"]) > 0
        # The Lean theorem's seed hypotheses at the same (h, β₀) follow from the certificate.
        point = only(filter(p -> p.h == h, lean))
        @test point.β == β0
        @test parse_rational(cert["seed"]["eps_s"]) <= point.s
        @test parse_rational(cert["seed"]["eps_t"]) <= point.t
    end
end

@testset "provenance and summary" begin
    for h in FIELDS
        meta = read_exact_json(certificate_path(h))["metadata"]
        @test verify_model_metadata(meta, ROOT; SOURCES...) == 3
    end
    summary = read_exact_json(joinpath(ROOT, "results", "summary.json"))
    @test [parse_rational(r["h"]) for r in summary["points"]] == FIELDS
    for r in summary["points"]
        c = read_exact_json(certificate_path(parse_rational(r["h"])))
        @test all(r[k] == c["seed"][k] for k in ("n0", "beta0", "eps_s", "eps_t"))
        @test r["gap_lower_bound"] == c["gap_lower_bound"]
    end
end
