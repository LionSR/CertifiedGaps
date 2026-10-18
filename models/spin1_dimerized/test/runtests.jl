# Fast tests: the exact frame change from the physical chain to the rotated-frame RingModel, the
# exact model, reflection and RP data, the core spectrum against an independent dense
# physical-frame Hamiltonian, and every stored certificate (δ = 9/10 including a fresh core
# replay of its spectra), with tampering rejected.
using Test, LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate, CertifiedGaps.Results, CertifiedGaps.RingModels
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum, float_ring_spectrum
using CertifiedGaps.ExactArithmetic: parse_rational
using CertifiedGaps.Spins: integer_cartesian_generators
const D = DimerizedCertificate

"Exact operator A ⊗ A on sites `i` and `k` of `n` spin-1 sites (no matrix products)."
pair_op(A, i, k, n) = foldl(kron, [j in (i, k) ? Rational{Int}.(A) : Matrix{Rational{Int}}(I, 3, 3)
                                   for j in 1:n])

"Physical H on 2·cells spins in the S^z basis, written site by site: bond (j, j+1) has
coupling 1 + (−1)^j δ. Independent of the `RingModel` cell construction."
function physical_ring(δ, cells)
    n = 2cells
    H = zeros(Rational{Int}, 3^n, 3^n)
    for j in 1:n, (c, A) in ((1 // 2, D.XP), (-1 // 2, D.AP), (1, D.SZ))
        H += (1 + (-1)^j * δ) * c * pair_op(A, j, mod1(j + 1, n), n)
    end
    H
end

@testset "Exact frame unitary" begin
    V, Sx, A, Sz = D.V, D.XP, D.AP, D.SZ
    @test V * V' == I && V == V'
    # π rotation about y: S^x, S^z reverse, iS^y is kept.
    @test V * Sx * V' == -Sx && V * A * V' == A && V * Sz * V' == -Sz
    # Spin algebra in this basis: S·S has eigenvalues −2, −1, 1 (total spin 0, 1, 2).
    @test sort(round.(Int, eigvals(Float64.(D.physical_bond(1))))) == [-2; fill(-1, 3); fill(1, 5)]
    for J in (1 // 10, 19 // 10)
        @test kron(D.ONE, V) * D.physical_bond(J) * kron(D.ONE, V)' == D.rotated_bond(J)
        @test kron(V, D.ONE) * D.physical_bond(J) * kron(V, D.ONE)' == D.rotated_bond(J)
    end
    for δ in D.DELTAS
        @test is_frame_change(D.physical(δ), D.dimerized(δ), D.cell_unitary())
    end
    # The cell model of the physical chain is the site-by-site chain (ring of two cells).
    for δ in (9 // 10, 7 // 10)
        m = D.physical(δ)
        o, h = Rational{Int}.(real(m.onsite)), Rational{Int}.(real(m.bond))
        one9 = Matrix{Rational{Int}}(I, 9, 9)
        S = swap_matrix(9)
        @test kron(o, one9) + kron(one9, o) + h + S * h * S == physical_ring(δ, 2)
    end
end

# The Cartesian-frame identities behind METHOD.md's RP and X_β = K K* argument.
@testset "Cartesian frame identities" begin
    G = integer_cartesian_generators()
    S = [-im * g for g in G]                                          # S^α = −i G_α
    @test S[1] * S[2] - S[2] * S[1] == im * S[3]
    @test sum(s * s for s in S) == 2I
    @test -sum(kron(s, s) for s in S) == sum(kron(g, g) for g in G)
    F, one3 = swap_matrix(3), Matrix{Int}(I, 3, 3)
    @test all(F * kron(one3, g) * F == kron(g, one3) for g in G)     # Θ(A_α) = B_α
end

@testset "Exact model, reflection and RP data" begin
    for δ in D.DELTAS
        m = D.dimerized(δ)                        # charge, flip, reflection checked exactly
        @test m.flip == kron(D.P, D.P) && m.charge == [0, 1, 2, -1, 0, 1, -2, -1, 0]
        @test m.reflection == D.SWAP * kron(D.P, D.P)
        cert = D.rp(m, δ)
        @test verify_rp(m, cert) && encode_rp(cert)["form"] == "reflected"
    end
    m = D.dimerized(9 // 10)
    # The bare cell swap is not a lattice reflection here: it reverses the charge.
    @test_throws ArgumentError RingModel(9, m.bond, m.onsite; charge = m.charge,
        flip = m.flip, reflection = D.SWAP)
    @test verify_rp(m, rp_certificate(m; reflection = D.SWAP))       # exact LDL† with R
    @test_throws ArgumentError rp_certificate(m)                      # not RP without R
end

@testset "Core spectrum against the physical frame" begin
    for δ in (9 // 10, 7 // 10)
        reference = eigvals(Symmetric(Float64.(physical_ring(δ, 3))))
        spec = verified_ring_spectrum(D.dimerized(δ), 3)
        levels = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spec]))
        @test maximum(abs, levels - reference) < 1e-12
        @test maximum(abs, float_ring_spectrum(D.dimerized(δ), 3).E - reference) < 1e-12
        @test all(s.label["charge"] >= 0 for s in spec)             # ±q paired by the flip
    end
end

@testset "Stored certificates" begin
    paths = D.result_paths()
    @test basename.(paths) == ["certificate_d7.json", "certificate_d9.json"]
    for path in paths
        @test D.check(path)["gap_lower_bound"] == read_result(path)["gap_lower_bound"]
    end
    # δ = 9/10: a fresh core replay of the spectra of 2 and 4 cells (δ = 7/10: check.jl --full).
    path = D.certificate_path(9 // 10)
    cert = read_result(path)
    r = D.check(path; full = true)
    @test r["defect_s_upper"] <= parse_rational(cert["seed"]["eps_s"])
    @test r["defect_t_upper"] <= parse_rational(cert["seed"]["eps_t"])
    @test D.lean_seed(9 // 10) !== nothing && D.lean_seed(7 // 10) === nothing
    @test read_result(D.certificate_path(7 // 10))["seed"]["positive_transfer"] ==
          D.POSITIVE_TRANSFER
end

@testset "Tampering is rejected" begin
    path = D.certificate_path(9 // 10)
    original = read_result(path)
    function rejects(edit!)
        cert = deepcopy(original)
        edit!(cert)
        tmp = joinpath(mktempdir(), basename(path))
        write_result(tmp, cert)
        try
            D.check(tmp)
            return false
        catch err
            return err isa ErrorException
        end
    end
    halve(x) = string(parse_rational(x) / 2)
    @test rejects(c -> c["seed"]["eps_s"] = halve(c["seed"]["eps_s"]))
    @test rejects(c -> c["seed"]["beta0"] = "4")
    @test rejects(c -> c["gap_lower_bound"] = halve(c["gap_lower_bound"]))
    @test rejects(c -> c["params"]["delta"] = "4//5")
    @test rejects(c -> c["interval_bound_decimal"] = 2 * c["interval_bound_decimal"])
    @test rejects(c -> delete!(c["model"], "reflection"))
    @test rejects(c -> (h = c["metadata"]["source_sha256"]; h[first(keys(h))] = repeat("0", 64)))
end
