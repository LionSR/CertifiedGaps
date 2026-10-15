# Checker for the dimerised Heisenberg certificates (about 10 s). From the stored files it
#   1. re-derives every exact RP identity of results/rp_audit.json at all six deltas;
#   2. for each results/verified_delta_*.json.gz, recomputes every block's scalar spectral
#      witness (residual, orthogonality, eigenvalue radius), the three outward log-partition
#      enclosures, the exact rational seeds and `certified_gap_bound`, and compares them exactly;
#   3. confirms that results/certificate.json copies those seeds exactly;
#   4. checks the recorded source hashes (src/provenance.jl).
# The fast check does not redo the eigenvector solves. `check.jl --full` (about 2.5 min) also
# re-solves every block and proves the stored seeds again (check_fresh_replay);
# scripts/independent_review.jl adds the stricter interval-containment test.
using CertifiedGaps, CertifiedGaps.Results, JSON3, Test, IntervalArithmetic, LinearAlgebra
isdefined(@__MODULE__, :check_historical_hashes) ||
    include(joinpath(@__DIR__, "src", "provenance.jl"))
isdefined(@__MODULE__, :VerifiedSpectra) ||
    include(joinpath(@__DIR__, "src", "verified_spectra.jl"))
isdefined(@__MODULE__, :DimerizedRP) || include(joinpath(@__DIR__, "src", "rp.jl"))
using .VerifiedSpectra, .DimerizedRP
BLAS.set_num_threads(1)

const RESULTS_DIR = joinpath(@__DIR__, "results")

"Parse an exact `\"num//den\"` string."
parseq(s) = (p = split(String(s), "//"); parse(BigInt, p[1]) // parse(BigInt, p[2]))

"Re-derive the exact RP audit (reflected form holds, identical-factor form fails)."
function check_audit()
    data = read_exact_json(joinpath(RESULTS_DIR, "rp_audit.json"))
    @test data["exact_algebraic_audit"] === true
    # The audit is about the RP hypothesis only; the seeds live in the verified certificates.
    @test data["gap_certificate"] === false
    @test data["metadata"]["certified"] === false
    @test length(data["parameters"]) == 6
    for row in data["parameters"]
        a = audit_rp(parseq(row["delta"]))
        @test a.defect_norm2 == parseq(row["swap_defect_frobenius_squared"])
        @test a.real_rhs_error_norm2 == parseq(row["proposed_real_rhs_error_frobenius_squared"])
        @test a.defect_norm2 > 0
        @test row["conjugate_identity_exact"] === true
        @test a.reflected_real_identity === true
        @test row["reflected_real_identity_exact"] === true
        @test row["theta_invariant_onsite_exact"] === true
        @test row["repository_general_RP_eligible"] === true
        @test row["real_identical_factor_identity"] === false
        @test row["rational_seeds"] === nothing
        @test row["certified_gap_bound"] === nothing
    end
    check_source_hashes(data["metadata"])
end

"Recompute every witness, partition enclosure, seed and the bound of one verified certificate."
function check_verified_file(path)
    data = JSON3.read(read_json_text(path))
    @test data.original_stronger_RP_condition_satisfied === false
    @test data.general_reflected_RP_exact === true
    @test 0 <= parseq(data.delta) < 1
    @test audit_rp(parseq(data.delta)).reflected_real_identity
    check_source_hashes(data.metadata)
    spectra = Dict{Int,Tuple{Vector{Float64},Vector{Rational{BigInt}}}}()
    for ring in data.rings
        E = Float64[]; radii = Rational{BigInt}[]; totaldim = 0
        for block in ring.blocks
            d = Int(block.dimension)
            w = VerifiedSpectra.scalar_witness(d, Int(ring.spins), parseq(block.h_component_max),
                parseq(block.q_component_max), parseq(block.lambda_max_abs),
                parseq(block.raw_residual_component_max),
                parseq(block.raw_orthogonality_component_max), parseq(block.matrix_error_upper))
            @test w.steps * VerifiedSpectra.UNIT < 1
            @test w.ebound == parseq(block.orthogonality_norm_upper) < 1
            @test w.rbound == parseq(block.residual_norm_upper)
            @test w.radius == parseq(block.eigenvalue_radius) > 0
            vals = parseq.(block.eigenvalues)
            @test length(vals) == d
            @test all(x -> Rational{BigInt}(Float64(x)) == x, vals)
            @test issorted(vals)
            for _ in 1:Int(block.multiplicity)
                append!(E, Float64.(vals)); append!(radii, fill(w.radius, d))
            end
            totaldim += d * Int(block.multiplicity)
        end
        # Completeness: every state of the ring is counted once, sector by sector.
        @test totaldim == 2^Int(ring.spins) == Int(ring.dimension) == length(E)
        for q in 0:Int(ring.spins)÷2
            @test sum(Int(b.dimension) for b in ring.blocks if Int(b.down_spins) == q) ==
                binomial(Int(ring.spins), q)
        end
        spectra[Int(ring.spins)] = (E, radii)
    end
    beta = parseq(data.beta0)
    small = enclosed_logZ(spectra[8]..., beta)
    large = enclosed_logZ(spectra[16]..., beta)
    double = enclosed_logZ(spectra[16]..., 2beta)
    for (key, value) in ((:n_beta, small), (:var"2n_beta", large), (:var"2n_2beta", double))
        stored = parseq.(data.log_partition_enclosures[key])
        @test stored == [rational_lower(value), rational_upper(value)]
    end
    es = rational_upper(1 - exp(large - 2small)); et = rational_upper(1 - exp(double - 2large))
    @test es == parseq(data.epsilon_s_upper)
    @test et == parseq(data.epsilon_t_upper)
    result = certified_gap_bound(es, et, beta)
    @test result.bound !== nothing
    @test Rational{BigInt}(result.bound) == parseq(data.recurrence_bound_exact_binary)
    @test result.bound == data.recurrence_bound
    println("verified stored delta=$(data.delta) bound=$(result.bound)")
end

"The compressed verified certificates, sorted by file name."
verified_paths() = sort(filter(p -> startswith(basename(p), "verified_delta_") &&
    endswith(p, ".json.gz"), readdir(RESULTS_DIR; join=true)))

"Check all three verified certificates and the plain summary."
function check_verified_results()
    paths = verified_paths()
    @test length(paths) == 3
    for path in paths
        check_verified_file(path)
    end
    check_certificate_summary(paths)
end

"The plain results/certificate.json must copy the seeds of the compressed certificates exactly."
function check_certificate_summary(paths)
    summary = JSON3.read(read(joinpath(RESULTS_DIR, "certificate.json"), String))
    check_source_hashes(summary.metadata)
    rows = Dict(String(r.source_file) => r for r in summary.certificates)
    @test length(rows) == length(paths)
    for path in paths
        r = rows["results/" * basename(path)]; data = JSON3.read(read_json_text(path))
        @test String(r.source_file_sha256) == sha256_file(path)
        @test Int(r.n0_cells) == Int(data.n0_cells) == 4 && iseven(Int(r.n0_cells))
        for key in (:delta, :beta0, :epsilon_s_upper, :epsilon_t_upper,
                :recurrence_bound_exact_binary)
            @test String(r[key]) == String(data[key])
        end
        for (k, v) in pairs(data.log_partition_enclosures)
            @test String.(r.log_partition_enclosures[k]) == String.(v)
        end
    end
end

"""
    refined_witness(rec, L)

Rigorous residual, orthogonality and eigenvalue-radius bounds for one freshly solved block
record `rec` of an `L`-spin ring, sharper than `VerifiedSpectra.scalar_witness` (Cauchy on the
concatenated real dot products). Used by `check_fresh_replay` and scripts/independent_review.jl.
"""
function refined_witness(rec, L)
    d = Int(rec["dimension"])
    h = parseq(rec["matrix_error_upper"])
    q = parseq(rec["q_component_max"])
    lam = parseq(rec["lambda_max_abs"])
    rawr = parseq(rec["raw_residual_component_max"])
    rawe = parseq(rec["raw_orthogonality_component_max"])
    coarse = VerifiedSpectra.scalar_witness(d, L, parseq(rec["h_component_max"]),
        q, lam, rawr, rawe, h)
    @assert coarse.ebound < 1 # Thus ||Q|| < 2, independently of the refinement.
    a = big(3) * L // 4
    # Cauchy on two concatenated real dot-product vectors gives these sums.
    rerr = min(coarse.rerr, coarse.gamma * (2 * (a + h) + q * lam) + coarse.underflow)
    eerr = min(coarse.eerr, coarse.gamma * 5 + coarse.underflow)
    # ||complex dxd matrix||F <= sqrt(2)*d*componentmax < (3/2)*d*componentmax.
    # This also strictly improves the coarse radius of scalar/small blocks.
    factor = big(3) * d // 2
    r = factor * (rawr + rerr); e = factor * (rawe + eerr)
    @assert e < 1
    return (; r, e, radius = r + 2 * h + (a + lam) * e)
end

"""
    check_fresh_replay()

Slow replay (`check.jl --full`, about 2.5 min): re-solve both rings of every certified delta
from the model, enclose every eigenvalue with `refined_witness`, and require the fresh rational
seeds to be at most the stored ones, so the fresh data prove the stored bound again. Unlike
scripts/independent_review.jl it does not require each fresh interval to lie inside the stored
one; that stricter test can fail on another BLAS without affecting validity (README).
"""
function check_fresh_replay()
    for path in verified_paths()
        data = JSON3.read(read_json_text(path))
        delta = parseq(data.delta); beta = parseq(data.beta0)
        fresh = Dict{Int,Tuple{Vector{Float64},Vector{Rational{BigInt}}}}()
        for stored in data.rings
            L = Int(stored.spins)
            ring, _, _ = verified_ring(delta, L)
            @test ring["dimension"] == stored.dimension
            @test length(ring["blocks"]) == length(stored.blocks)
            centers = Float64[]; radii = Rational{BigInt}[]
            for (a, b) in zip(ring["blocks"], stored.blocks)
                for key in ("dimension", "down_spins", "cell_momentum", "spin_flip_parity",
                        "multiplicity")
                    @test a[key] == b[Symbol(key)]
                end
                w = refined_witness(a, L)
                vals = Float64.(parseq.(a["eigenvalues"]))
                for _ in 1:Int(a["multiplicity"])
                    append!(centers, vals); append!(radii, fill(w.radius, length(vals)))
                end
            end
            @test length(centers) == 2^L
            fresh[L] = (centers, radii)
        end
        small = enclosed_logZ(fresh[8]..., beta); large = enclosed_logZ(fresh[16]..., beta)
        double = enclosed_logZ(fresh[16]..., 2beta)
        es = rational_upper(1 - exp(large - 2small)); et = rational_upper(1 - exp(double - 2large))
        @test es <= parseq(data.epsilon_s_upper)
        @test et <= parseq(data.epsilon_t_upper)
        bound = certified_gap_bound(es, et, beta).bound
        @test bound !== nothing && bound >= data.recurrence_bound
        println("fresh replay delta=$(data.delta) stored=$(data.recurrence_bound) fresh=$bound")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Dimerised Heisenberg certificate" begin
        @testset "Exact RP audit" begin
            check_audit()
        end
        @testset "Verified spectra, partition sums, seeds and bounds" begin
            check_verified_results()
        end
        if "--full" in ARGS
            @testset "Fresh spectral replay of all three seeds" begin
                check_fresh_replay()
            end
        end
    end
end
