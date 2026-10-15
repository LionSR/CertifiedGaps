# Fast tests (about 2 min): the fast checker (check.jl) plus independent validations of the RP
# factors, the sector code, the enclosures on small rings, the conjugate-momentum reuse, the
# dimer perturbation bound, the recurrence rejection and the checker's rejection of corrupted
# certificates.
include(joinpath(@__DIR__, "..", "check.jl"))
include(joinpath(@__DIR__, "..", "src", "reflection_conventions.jl"))
using CertifiedGaps.RingModels
const V = Certification.VerifiedSpectra
const Q = Rational{BigInt}

@testset "Exact RP audit (the RP hypothesis; seeds are checked below)" begin
    check_audit()
end

@testset "Exact alternative cell conventions" begin
    conventions = ReflectionConventions.audit_conventions()
    @test length(conventions) == 2
    @test all(row.reflected_form && !row.uniform_identical_form for row in conventions)
end

@testset "Exact conjugated-factor RP and palindromic cell" begin
    for delta in DELTAS
        a = audit_rp(delta)
        model = RingModel(4, a.transformed_bond, a.transformed_onsite)
        rp = rp_certificate(model; factors=a.channels, weights=fill(1 - delta, 3))
        @test verify_rp(model, rp)
        @test isdiag(a.transformed_onsite)
        @test all(iszero, imag.(a.transformed_onsite))
        @test a.basis_numerator' * a.basis_numerator == 2Matrix{Complex{Q}}(I, 4, 4)
        wrong = RingModel(4, -a.transformed_bond, a.transformed_onsite)
        @test !verify_rp(wrong, rp)
    end
end

include("verified_independent.jl")

@testset "Exact rational residuals of every 8-spin block" begin
    CQ = Complex{Q}
    for delta in DELTAS[1:3]
        blocks = 0
        for q in 0:4
            orbits, flip = V.group_orbits(8, q)
            for k in 0:3, p in (flip ? (1, -1) : (1,))
                built = V.verified_block(delta, 8, q, k, p, orbits)
                built === nothing && continue
                H, err = built
                record, _, _ = V.enclose_block(H, err, 8)
                dec = eigen(Hermitian(H))
                U = CQ.(dec.vectors)
                R = CQ.(H) * U - U * Diagonal(Q.(dec.values))
                G = U' * U - I
                @test sum(abs2, R) <= parseq(record["residual_norm_upper"])^2
                @test sum(abs2, G) <= parseq(record["orthogonality_norm_upper"])^2
                blocks += 1
            end
        end
        @test blocks > 0
    end
end

@testset "Scalable dimer perturbation and recurrence rejection" begin
    d = dimer_bounds(999999 // big(1000000), 128, 8 // big(1))
    # A genuine beyond-ED example: largest ring 512 spins, no diagonalisation.
    b = certified_gap_bound(Q(sup(d.s)), Q(sup(d.t)), 8 // big(1)).bound
    @test b !== nothing && b > 0
    for delta in DELTAS
        bad = dimer_bounds(delta, 4, 6 // big(1))
        @test inf(bad.s) > 1 // 2
        @test Certification.rejects(Q(inf(bad.s)), Q(inf(bad.t))) !== nothing
    end
    @test Certification.rejects(1 // big(2), 0 // big(1)) == 0
    @test Certification.rejects(1 // big(100), 1 // big(100)) === nothing
    @test_throws ErrorException seed_intervals(Dict(), 5, 6 // big(1))
end

@testset "Exact conjugate momentum reuse" begin
    # verified_ring reuses the spectrum of (k, p) for (n-k, p); the blocks must be conjugate.
    for (L, charges) in ((8, 0:4), (16, (2, 8))), q in charges
        n = L ÷ 2
        orbits, flip = V.group_orbits(L, q)
        for k in 1:(n - 1)÷2, p in (flip ? (1, -1) : (1,))
            left = V.verified_block(3 // big(10), L, q, k, p, orbits)
            right = V.verified_block(3 // big(10), L, q, n - k, p, orbits)
            @test (left === nothing) == (right === nothing)
            left === nothing && continue
            @test left[1] == conj.(right[1])
            @test left[2] == right[2]
        end
    end
end

include("numerical_blocks.jl")
include(joinpath(@__DIR__, "..", "src", "Diagnostics.jl"))
@testset "Sparse momentum/parity blocks against complete spectra" begin
    for L in (4, 6, 8, 10)
        reference, _ = dense_spectrum(0.3, L)
        sparse, records = Diagnostics.low_spectrum(0.3, L; dense_max=2000)
        @test sparse.exact
        @test maximum(abs, sparse.E - reference.E) < 1e-11
        @test sum(parse(Int, b["dimension"]) * parse(Int, b["multiplicity"])
            for b in records) == 2^L
    end
end

@testset "Checker rejects corrupted theorem inputs" begin
    record = load_record(seed_path(1 // big(2)))
    spectra = stored_spectra(record)
    @test check_record(record; spectra, check_hashes=false) === spectra
    for (key, value) in (("eps_s", "0//1"), ("eps_t", "0//1"),
        ("gap_lower", "1//1"), ("n0_cells", "5"), ("validity_sites", "all even N >= 8"))
        corrupted = deepcopy(record)
        corrupted["seeds"][1][key] = value
        @test_throws ErrorException check_record(corrupted; spectra, check_hashes=false)
    end
    corrupted = deepcopy(record)
    empty!(corrupted["metadata"]["source_sha256"])
    @test_throws ErrorException check_record(corrupted; spectra)
    corrupted = deepcopy(record)
    popfirst!(corrupted["negative_results"][1]["cover"])
    @test_throws ErrorException check_record(corrupted; spectra, check_hashes=false)
    corrupted = deepcopy(record)
    empty!(corrupted["negative_results"])
    @test_throws ErrorException check_record(corrupted; spectra, check_hashes=false)
    corrupted = deepcopy(record)
    corrupted["reflection_positivity"]["factors"][1]["weight"] = "-1//1"
    @test_throws ErrorException check_record(corrupted; spectra, check_hashes=false)
    corrupted = deepcopy(record)
    corrupted["dimer_perturbation"]["norm_small"] = "0//1"
    @test_throws ErrorException check_record(corrupted; spectra, check_hashes=false)
    # Stored witnesses: a smaller radius or residual maximum must be caught by the fast check.
    for key in ("eigenvalue_radius", "raw_residual_component_max")
        corrupted = deepcopy(record)
        block = argmax(b -> Certification.asint(b["dimension"]), corrupted["rings"][end]["blocks"])
        block[key] = string(exact(block[key]) / 2)
        @test_throws ErrorException stored_spectra(corrupted)
    end
    corrupted = deepcopy(record)
    pop!(corrupted["rings"][end]["blocks"])
    @test_throws ErrorException stored_spectra(corrupted)
end

@testset "Seed certificates and exclusions (stored spectra)" begin
    check_seed_certificates()
end

@testset "Summary index" begin
    check_certificate_summary()
end
