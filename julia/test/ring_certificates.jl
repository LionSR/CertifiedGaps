# Generic ring models, verified spectra and end-to-end gap certificates.

using CertifiedGaps.RingModels, CertifiedGaps.VerifiedSpectra, CertifiedGaps.GapCertificates
const JSON3 = CertifiedGaps.Results.JSON3          # not a test dependency of its own

const HALF = 1 // 2
const SX, ISY, SZ = [0 HALF; HALF 0], [0 HALF; -HALF 0], [HALF 0; 0 -HALF]

# Staggered-field Heisenberg after the sublattice rotation: −SˣSˣ + SʸSʸ − SᶻSᶻ − h_s Sᶻ.
rotated_staggered(hs) = RingModel(2, -kron(SX, SX) - kron(ISY, ISY) - kron(SZ, SZ), -hs * SZ;
    charge = [0, 1], modulus = 2, name = "staggered Heisenberg, rotated frame")

levels(spec) = sort(reduce(vcat, [repeat(s.centers, s.multiplicity) for s in spec]))

@testset "ring models and exact RP certificates" begin
    m = rotated_staggered(1)
    cert = rp_certificate(m; factors = [SX, ISY, SZ])
    @test verify_rp(m, cert) && cert.weights == [1, 1, 1]
    @test encode_rp(cert)["factors"][2]["transpose_symmetry"] == "antisymmetric"
    @test encode_rp(cert)["form"] == "real"
    # The exact LDL† factorisation finds an equivalent certificate on its own.
    @test verify_rp(m, rp_certificate(m))
    # Ferromagnetic-sign Heisenberg in the physical frame is not RP in the standard basis.
    @test_throws ArgumentError rp_certificate(RingModel(2, kron(SX, SX) + kron(SZ, SZ)))
    @test_throws ArgumentError RingModel(2, -kron(SX, SX); charge = [1, -1])
    # Exact JSON round trip, through the parser that reads integers as strings.
    record = parse_exact_json(JSON3.write(Dict("m" => encode_model(m), "rp" => encode_rp(cert))))
    @test decode_model(record["m"]).bond == m.bond
    @test verify_rp(decode_model(record["m"]), decode_rp(record["rp"]))
end

@testset "reflected RP form: dimerised Heisenberg in two-site cells" begin
    # Cell (a, b): rotated-frame bond J₁ inside the cell (on-site term), J₂ from b to the next a.
    rotated(J) = J * (-kron(SX, SX) - kron(ISY, ISY) - kron(SZ, SZ))
    one = [1 0; 0 1]
    J1, J2 = 3 // 2, 1 // 2
    bond = kron(one, rotated(J2), one)          # site b of the left cell, site a of the right
    swap = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    cell = RingModel(4, bond, rotated(J1); name = "dimerised Heisenberg cells")
    # −h = J₂ Σ (1 ⊗ S) ⊗ (S ⊗ 1), and S ⊗ 1 = R (1 ⊗ S) R with R the swap of the cell.
    factors = [kron(one, A) for A in (SX, ISY, SZ)]
    cert = rp_certificate(cell; factors, weights = fill(J2, 3), reflection = swap)
    @test verify_rp(cell, cert)
    record = encode_rp(cert)
    @test record["form"] == "reflected" && decode_rp(record).reflection == swap
    @test verify_rp(cell, rp_certificate(cell; reflection = swap))     # exact LDL† with R
    # Without the reflection the bond is not RP in this basis.
    @test_throws ArgumentError rp_certificate(cell)
    # A reflection that does not commute with the on-site term is refused.
    @test !verify_rp(RingModel(4, bond, kron(SZ, one)), RPCertificate(cert.weights,
        cert.factors, cert.reflection))
    # Spectrum against ED of the alternating chain in the physical frame.
    heis = real(kron(SZ, SZ) + kron(SX, SX) - kron(ISY, ISY))
    chain = Chain(2, L -> [(j, mod1(j + 1, L), Float64.(isodd(j) ? J1 : J2) * heis) for j in 1:L])
    E = sort(spectrum(chain, 8; dense_max = 10^6).E)
    @test maximum(abs, levels(verified_ring_spectrum(cell, 4)) - E) < 1e-12
end

@testset "verified spectra" begin
    # Against independent dense ED in the physical frame (unitarily equivalent for even L).
    m = rotated_staggered(1)
    for L in (4, 6)
        spec = verified_ring_spectrum(m, L)
        E = sort(spectrum(staggered_heisenberg(1).chain, L; dense_max = 10^6).E)
        @test maximum(abs, levels(spec) - E) < 1e-12
        @test maximum(s.radius for s in spec) < 1e-11
    end
    # Transverse-field Ising: log Z enclosure contains the exact free-fermion value.
    t = RingModel(2, -kron(2SZ, 2SZ), -2 * (2SX))
    z = log_partition_enclosure(verified_ring_spectrum(t, 10), 3 // 2)
    @test in_interval(Float64(tfim_logZ(10, 3 // 2, 2)), z) && diam(z) < 1e-9
    # A complex model without charge, all momenta built separately.
    B = [1 0 0 0; 0 0 (1 + 2im)//3 0; 0 (1 - 2im)//3 0 0; 0 0 0 -HALF]
    O = [0 HALF*im; -HALF*im 0]
    c = RingModel(2, B, O)
    E = sort(spectrum(Chain(2, ring_bonds(ComplexF64.(B)); onsite = ComplexF64.(O)), 5;
        dense_max = 10^6).E)
    @test maximum(abs, levels(verified_ring_spectrum(c, 5)) - E) < 1e-12
    # The enclosure accounts for a perturbation of the matrix.
    M = [2.0 1.0; 1.0 2.0]
    centers, r = enclose_hermitian(M; matrix_error = 1 // 1000)
    @test centers ≈ [1.0, 3.0] && 1 // 1000 < r < 1 // 1000 + 1 // 10^12
    # Tridiagonal Toeplitz blocks with exact spectrum cos(kπ/(n+1)); the complex block is
    # diagonally unitarily equivalent (off-diagonal i/2). The radius must cover the truth.
    n = 300
    exact = sort([cospi(big(k) / (n + 1)) for k in 1:n])
    for off in (0.5, 0.5im)
        T = diagm(1 => fill(off, n - 1), -1 => fill(conj(off), n - 1))
        e = verified_eigen(off isa Complex ? T : real(T))
        @test maximum(abs.(exact .- e.values)) <= e.radius < 1 // 10^10
        # The polar-factor radius is kept: it halves the orthogonality term of Ostrowski's.
        η, r, λ = e.orthogonality, e.residual, maximum(abs, e.values)
        @test e.radius < (2η * λ + (1 + η) * r) / (1 - η)
    end
    @test_throws ArgumentError verified_eigen([1.0 2.0; nextfloat(2.0) 1.0])
end

@testset "charge flip: sectors q and −q built once" begin
    # Physical-frame Heisenberg (SʸSʸ = −(iSʸ)(iSʸ)) with the U(1) charge 2Sᶻ, reversed by 2Sˣ.
    bond = kron(SX, SX) - kron(ISY, ISY) + kron(SZ, SZ)
    plain = RingModel(2, bond; charge = [1, -1])
    flipped = RingModel(2, bond; charge = [1, -1], flip = 2SX)
    for L in (6, 7)
        a, b = verified_ring_spectrum(plain, L), verified_ring_spectrum(flipped, L)
        @test maximum(abs, levels(a) - levels(b)) < 1e-12
        @test length(b) < length(a) && all(s.label["charge"] >= 0 for s in b)
        @test inf(log_partition_enclosure(a, 2)) <= sup(log_partition_enclosure(b, 2)) &&
              inf(log_partition_enclosure(b, 2)) <= sup(log_partition_enclosure(a, 2))
    end
    # A field breaks the flip; a flip that does not reverse the charge is refused.
    @test_throws ArgumentError RingModel(2, bond, SZ; charge = [1, -1], flip = 2SX)
    @test_throws ArgumentError RingModel(2, bond; charge = [1, -1], flip = [1 0; 0 1])
    # The flip is part of the exact model record and survives a certificate round trip.
    m = RingModel(2, bond - kron(SZ, SZ) / 2; charge = [1, -1], flip = 2SX)
    @test decode_model(parse_exact_json(JSON3.write(encode_model(m)))).flip == m.flip
end

@testset "end-to-end gap certificate" begin
    # A fast seed deep in the admissible region: field h_s = 6, (n₀, β₀) = (4, 1).
    m = rotated_staggered(6)
    cert = certify_gap(m; n0 = 4, β0 = 1, rp = rp_certificate(m; factors = [SX, ISY, SZ]))
    bound = parse_rational(cert["gap_lower_bound"])
    @test 0 < bound < cert["interval_bound_decimal"]
    # The exact chain is not capped at 1 and stays close to the interval bound.
    @test 1 < bound && cert["interval_bound_decimal"] - bound < 1e-2
    stored = parse_exact_json(JSON3.write(cert))
    @test check_certificate(stored)["gap_lower_bound"] == cert["gap_lower_bound"]
    # A tampered seed defect is rejected.
    stored["seed"]["eps_s"] = string(parse_rational(stored["seed"]["eps_s"]) / 2)
    @test_throws ErrorException check_certificate(stored)
    # An inadmissible seed is refused.
    @test_throws ErrorException certify_gap(rotated_staggered(1); n0 = 2, β0 = 1)
end

@testset "blocking, seed reuse, odd seeds, gaps and stored spectra" begin
    # The rotated staggered chain conserves the sublattice charge z_a − z_b only in 2-site cells.
    m = rotated_staggered(6)
    rp = rp_certificate(m; factors = [SX, ISY, SZ])
    blocking = (; cell = 2, charge = [0, 1, -1, 0])
    cells = cell_model(m, 2; blocking.charge)
    @test levels(verified_ring_spectrum(cells, 4)) ≈ levels(verified_ring_spectrum(m, 8))
    @test_throws ArgumentError cell_model(m, 2; charge = [0, 1, 1, 0])
    seed = seed_bounds(m, 4, 1; blocking)
    plain = certify_gap(m; n0 = 4, β0 = 1, rp)
    cert = certify_gap(m; n0 = 4, β0 = 1, rp, blocking, seed)
    @test cert["blocking"]["cell"] == 2 && cert["model"] == plain["model"]
    @test abs(parse_rational(cert["seed"]["eps_s"]) - parse_rational(plain["seed"]["eps_s"])) <
          1 // 10^9
    stored = parse_exact_json(JSON3.write(cert))
    @test check_certificate(stored)["gap_lower_bound"] == cert["gap_lower_bound"]
    # A seed for other arguments is refused, and so is a stored cell model of another chain.
    @test_throws ArgumentError certify_gap(m; n0 = 4, β0 = 2, rp, blocking, seed)
    stored["blocking"]["model"] = encode_model(cell_model(rotated_staggered(5), 2;
        blocking.charge))
    @test_throws ErrorException check_certificate(stored; spectra = false)
    # Reused spectra give the same seed.
    again = seed_bounds(m, 4, 1; blocking, spectra = Dict(4 => seed.small, 8 => seed.large))
    @test (again.eps_s, again.eps_t) == (seed.eps_s, seed.eps_t)
    # Odd seed lengths only with a positivity reference, which the certificate records.
    @test_throws ArgumentError seed_bounds(m, 3, 1)
    odd = certify_gap(m; n0 = 3, β0 = 1, rp, positive_transfer = "test: X_β ≥ 0 assumed")
    @test check_certificate(parse_exact_json(JSON3.write(odd)); spectra = false) isa Dict
    # Gap and ground-energy enclosures against dense ED of the 8-site ring.
    E = levels(seed.large)
    gl, gu = spectral_gap_enclosure(seed.large)
    el, eu = ground_energy_enclosure(seed.large)
    @test 0 < gl <= E[2] - E[1] <= gu && gu - gl < 1 // 10^9 && el <= E[1] <= eu
    @test spectral_gap_enclosure([BlockSpectrum([0.0, 1.0], 0, 2, Dict{String,Any}())]) == (0, 0)
    # Lossless storage, and log Z over a β interval.
    path = tempname() * ".json.gz"
    back = read_spectrum(write_spectrum(path, seed.large))
    @test [s.centers for s in back] == [s.centers for s in seed.large]
    @test [s.radius for s in back] == [s.radius for s in seed.large]
    box = log_partition_enclosure(back, interval(1.0, 1.5))
    @test issubset_interval(log_partition_enclosure(back, 5 // 4), box)
    s0 = minimum(E)
    @test isequal_interval(shifted_log_partition(back, 1, s0) - interval(s0),
        log_partition_enclosure(back, 1))
end
