# Fast unit tests of the model definitions in src/common.jl: the exact sublattice rotation to the
# RP frame, the RP data, the cell charge, and the verified spectra against an independent
# physical-frame exact diagonalization. The stored certificates are checked by check.jl.
using Test
include(joinpath(@__DIR__, "..", "src", "common.jl"))
using CertifiedGaps.ExactDiagonalization: spectrum
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum

@testset "exact sublattice rotation to the RP frame" begin
    @test U * U' == I2 && U * X * U' == -X && U * J * U' == J && U * Z * U' == -Z
    for h in FIELDS
        # Physical H = Σ S·S − h Σ_j (−1)^j Sᶻ_j in two-site cells (a, b) = (2x − 1, 2x): the
        # field is +hZ on a and −hZ on b. Rotating a by U gives the blocked rotated model on
        # every ring of cells, that is, on every even ring of sites.
        physical = RingModel(4, kron(I2, physical_bond(), I2),
            physical_bond() + h * kron(Z, I2) - h * kron(I2, Z))
        @test is_frame_change(physical, cell_model(site_model(h), 2), kron(U, I2))
    end
end

@testset "RP data, reflection and the two-site cell charge" begin
    for h in FIELDS
        m = site_model(h)
        @test verify_rp(m, site_rp(h)) && isreal_model(m) && m.reflection == I2
        # In cells (a, b) the charge Sᶻ_a − Sᶻ_b is minus the physical Sᶻ (U flips site a).
        @test BLOCKING.charge == [Int(Z[a, a] - Z[b, b]) for a in 1:2 for b in 1:2]
        @test cell_model(m, BLOCKING.cell; charge = BLOCKING.charge).d == 4
        # The cell swap mirrors the bond but reverses the cell charge, so the blocked model
        # cannot carry it as a reflection.
        swap = swap_matrix(2)
        @test_throws ArgumentError cell_model(m, 2; charge = BLOCKING.charge, reflection = swap)
        @test cell_model(m, 2; reflection = swap).reflection == swap
    end
end

@testset "verified spectra against physical-frame ED" begin
    ref(L) = spectrum(CertifiedGaps.Models.staggered_heisenberg(3 // 4).chain, L;
        dense_max = 10^6).E
    levels(spec) = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spec]))
    # The site model (Sᶻ parity, real reflection blocks) and the blocked cell model.
    for L in (4, 6)
        @test maximum(abs, levels(verified_ring_spectrum(site_model(3 // 4), L)) - ref(L)) <
              1e-12
    end
    cell = cell_model(site_model(3 // 4), BLOCKING.cell; charge = BLOCKING.charge)
    @test maximum(abs, levels(verified_ring_spectrum(cell, 4)) - ref(8)) < 1e-12
end
