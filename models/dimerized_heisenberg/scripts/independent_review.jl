# Independent acceptance: fresh solves, rigorous refined radii, stored-interval containment.
# Stricter than `check.jl --full`: every eigenvalue interval regenerated around the fresh centre
# with `refined_witness` must lie inside the stored interval. This does not require
# bitwise-identical eigenvectors, can fail on another BLAS without affecting validity (README),
# and never overwrites production data. About 4 minutes; writes results/acceptance.json
# (gitignored).
include(joinpath(@__DIR__, "..", "check.jl"))
const QQReview = Rational{BigInt}
reviewrecords = Any[]
elapsed = @elapsed @testset "Independent acceptance: containment from fresh full spectra" begin
    for delta in DELTAS
        path = seed_path(delta)
        record = read_result(path)
        println("ACCEPTANCE RERUN delta=", delta); flush(stdout)
        fresh = Dict{Int,Any}(); mismatches = Dict{String,Int}(); minmargin = nothing
        for stored in record["rings"]
            L = Certification.asint(stored["spins"])
            ring, _, _ = verified_ring(delta, L)
            @test length(ring["blocks"]) == length(stored["blocks"])
            centers = Float64[]; radii = QQReview[]
            for (a, b) in zip(ring["blocks"], stored["blocks"])
                for key in ("dimension", "down_spins", "cell_momentum", "spin_flip_parity",
                        "multiplicity")
                    @test string(a[key]) == string(b[key]) # stored integers may be strings
                end
                for key in ("raw_residual_component_max", "raw_orthogonality_component_max",
                        "eigenvalue_radius", "eigenvalues")
                    a[key] == b[key] || (mismatches[key] = get(mismatches, key, 0) + 1)
                end
                w = refined_witness(a, L)
                av = parse_rational.(a["eigenvalues"]); bv = parse_rational.(b["eigenvalues"])
                margin = parse_rational(b["eigenvalue_radius"]) - maximum(abs.(av .- bv)) - w.radius
                # Ordered Weyl intervals regenerated around the fresh centres lie wholly inside
                # the stored intervals, including multiplicities.
                @test margin >= 0
                minmargin = minmargin === nothing ? margin : min(minmargin, margin)
                for _ in 1:Int(a["multiplicity"])
                    append!(centers, Float64.(av)); append!(radii, fill(w.radius, length(av)))
                end
            end
            @test length(centers) == 2^L
            fresh[L] = (centers, radii)
        end
        bounds = Any[]
        seeds = copy(record["seeds"])
        record["survey_seed"] === nothing || push!(seeds, record["survey_seed"])
        for seed in seeds
            n = parse(Int, seed["n0_cells"]); beta = parse_rational(seed["beta0"])
            f = seed_intervals(fresh, n, beta)
            es = rational_upper(f.s); et = rational_upper(f.t)
            @test es <= parse_rational(seed["eps_s"]) && et <= parse_rational(seed["eps_t"])
            b = certified_gap_bound(es, et, beta).bound
            @test b !== nothing && b >= parse_rational(seed["gap_lower"])
            push!(bounds, Dict("n0_cells" => n, "beta0" => string(beta),
                "fresh_refined_bound" => b, "stored_bound" => seed["gap_lower"]))
        end
        push!(reviewrecords, Dict("delta" => string(delta), "certificate" => basename(path),
            "certificate_sha256" => sha256_file(path), "bounds" => bounds,
            "non_bitwise_equal_block_counts" => mismatches,
            "minimum_containment_margin" => string(minmargin)))
        println("ACCEPTED delta=", delta, " bounds=", bounds, " margin=", Float64(minmargin),
            " mismatch_counts=", mismatches); flush(stdout)
    end
end
root = normpath(joinpath(@__DIR__, "..", "..", ".."))
write_result(joinpath(@__DIR__, "..", "results", "acceptance.json"), Dict(
    "metadata" => run_metadata(Dict("scope" => "independent_containment_review");
        source_sha256=source_hashes(root; dirs=(), files=(
            "models/dimerized_heisenberg/scripts/independent_review.jl",
            "models/dimerized_heisenberg/check.jl",
            "models/dimerized_heisenberg/src/verified_spectra.jl")),
        packages=(CertifiedGaps, IntervalArithmetic), status="pass", certified=true),
    "backend_assumption" => "Classical binary64 CPU OpenBLAS GEMM with ordinary multiply/add/FMA",
    "seconds" => elapsed, "certificates" => reviewrecords))
