# The exact frame: Cartesian (Lean) basis to the Sᶻ basis of the core RingModel, the RP data, and
# the verified spectra against an independent dense Hamiltonian built in the Cartesian frame.
using Test, LinearAlgebra
using CertifiedGaps.RingModels: verify_rp, rp_certificate
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum
isdefined(@__MODULE__, :TransferCertificate) ||
    include(joinpath(dirname(@__DIR__), "src", "TransferCertificate.jl"))
using .TransferCertificate.LargeDModel

"Dense Σ S·S + D Σ (Sᶻ)² on the ring of L sites in the Cartesian basis (S^α)_{ab} = −i ε_{αab}."
function cartesian_ring(L, D)
    ε(a, b, c) = (a - b) * (b - c) * (c - a) ÷ 2
    S = [ComplexF64[-im * ε(α, a, b) for a in 1:3, b in 1:3] for α in 1:3]
    site(op, j) = kron(Matrix{ComplexF64}(I, 3^(j - 1), 3^(j - 1)), op,
        Matrix{ComplexF64}(I, 3^(L - j), 3^(L - j)))
    H = sum(site(s, j) * site(s, mod1(j + 1, L)) for j in 1:L for s in S)
    return H + Float64(D) * sum(site(S[3]^2, j) for j in 1:L)
end

@testset "Exact frame and reflection positivity" begin
    (; W, norms) = cartesian_frame()                  # asserts the exact identities
    @test W' * W == Diagonal(norms)
    for D in (4 // 1, 1603 // 500, 13 // 4)
        model = ring_model(D)
        @test verify_rp(model, rp_data(model))
        # The exact LDL factorization of the realignment finds the same reflected form.
        @test verify_rp(model, rp_certificate(model; reflection = REFLECTION))
    end
    # Without the reflection the bond is not reflection positive in the Sᶻ basis.
    @test_throws ArgumentError rp_certificate(ring_model(4 // 1))
end

@testset "Verified spectra against the Cartesian Hamiltonian, L ≤ 6" begin
    for L in 2:6, D in (4 // 1, 1603 // 500)
        spectrum = verified_ring_spectrum(ring_model(D), L)
        levels = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spectrum]))
        exact = eigvals(Hermitian(cartesian_ring(L, D)))
        radius = Float64(maximum(s.radius for s in spectrum))
        # The dense solve has its own rounding error of order 1e-13.
        @test length(levels) == 3^L && maximum(abs, levels - exact) <= radius + 1e-10
    end
end
