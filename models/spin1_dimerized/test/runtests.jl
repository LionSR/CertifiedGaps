# Fast tests: the exact frame unitary from the physical chain to the rotated-frame RingModel,
# the exact RP data, the floating-point scan decomposition, the core spectrum against the
# physical-frame scan code, and every stored certificate (δ = 9/10 including a fresh core
# replay of its spectra), with tampering rejected.
using Test, LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate, CertifiedGaps.Results, CertifiedGaps.RingModels
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum
using CertifiedGaps.ExactArithmetic: parse_rational
const D = DimerizedCertificate

"Exact operator A ⊗ A on sites `i` and `k` of `n` spin-1 sites (no matrix products)."
pair_op(A, i, k, n) = foldl(kron, [j in (i, k) ? Rational{Int}.(A) : Matrix{Rational{Int}}(I, 3, 3)
                                   for j in 1:n])

"Physical H on 2·cells spins in the S^z basis: bond (j, j+1) has coupling 1 + (−1)^j δ."
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
    # The whole ring of two cells: U H_phys U' equals the ring assembled from the RingModel.
    for δ in (9 // 10, 7 // 10)
        m = D.dimerized(δ)
        o, h = Rational{Int}.(real(m.onsite)), Rational{Int}.(real(m.bond))
        one9 = Matrix{Rational{Int}}(I, 9, 9)
        swap = zeros(Rational{Int}, 81, 81)
        for x in 1:9, y in 1:9
            swap[9(y - 1) + x, 9(x - 1) + y] = 1
        end
        H_model = kron(o, one9) + kron(one9, o) + h + swap * h * swap
        U = kron(D.cell_unitary(), D.cell_unitary())
        @test U * physical_ring(δ, 2) * U' == H_model
    end
end

@testset "Exact model and RP data" begin
    for δ in D.DELTAS
        m = D.dimerized(δ)                        # charge and flip conditions checked exactly
        @test m.flip == kron(D.P, D.P) && m.charge == [0, 1, 2, -1, 0, 1, -2, -1, 0]
        cert = D.rp(m, δ)
        @test verify_rp(m, cert) && encode_rp(cert)["form"] == "reflected"
    end
    m = D.dimerized(9 // 10)
    @test verify_rp(m, rp_certificate(m; reflection = D.SWAP))       # exact LDL† with R
    @test_throws ArgumentError rp_certificate(m)                      # not RP without R
end

@testset "Floating-point scan decomposition" begin
    delta = 1 // 2; cells = 2; groups = D.charge_groups(cells)
    reference = Float64[]
    for m in 0:2cells
        A, den = D.integer_hamiltonian(cells, delta, groups[m+1])
        E = eigvals(Symmetric(Matrix(A) / den))
        append!(reference, E); m > 0 && append!(reference, E)
    end
    spectrum, _, _ = D.dense_spectrum(cells, delta)
    @test length(spectrum.E) == 9^cells
    @test maximum(abs, spectrum.E - sort(reference)) < 1e-12
    # At the dimer limit each bond has energies −2, −1, 1 with degeneracy 1, 3, 5.
    decoupled, _, _ = D.dense_spectrum(cells, 1 // 1)
    single = [-4; fill(-2, 3); fill(2, 5)]
    @test maximum(abs, decoupled.E - sort([a + b for a in single for b in single])) < 1e-12
end

@testset "Core spectrum of the rotated frame against the physical frame" begin
    for δ in (9 // 10, 7 // 10)
        spec = verified_ring_spectrum(D.dimerized(δ), 3)
        levels = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spec]))
        reference, _, _ = D.dense_spectrum(3, δ)
        @test maximum(abs, levels - reference.E) < 1e-12
        @test all(s.label["charge"] >= 0 for s in spec)             # ±q paired by the flip
    end
end

@testset "Stored certificates" begin
    paths = D.result_paths()
    @test basename.(paths) == ["certificate_d7.json", "certificate_d9.json"]
    for path in paths
        @test D.check_point(path)["gap_lower_bound"] == read_result(path)["gap_lower_bound"]
    end
    # δ = 9/10: a fresh core replay of the spectra of 2 and 4 cells (δ = 7/10: check.jl --full).
    path = D.certificate_path(9 // 10)
    cert = read_result(path)
    r = D.check_point(path; full = true)
    @test r["defect_s_upper"] <= parse_rational(cert["seed"]["eps_s"])
    @test r["defect_t_upper"] <= parse_rational(cert["seed"]["eps_t"])
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
            D.check_point(tmp)
            return false
        catch err
            return err isa ErrorException
        end
    end
    halve(x) = string(parse_rational(x) / 2)
    @test rejects(c -> c["seed"]["eps_s"] = halve(c["seed"]["eps_s"]))
    @test rejects(c -> c["seed"]["beta0"] = "4")
    @test rejects(c -> c["gap_lower_bound"] = halve(c["gap_lower_bound"]))
    @test rejects(c -> c["delta"] = "4//5")
    @test rejects(c -> (h = c["metadata"]["source_sha256"]; h[first(keys(h))] = repeat("0", 64)))
end
