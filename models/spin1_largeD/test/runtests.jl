# Fast checker of the stored spin-1 large-D certificates (no dense diagonalisation).
#
# From the repository root:
#   OPENBLAS_NUM_THREADS=1 julia --project=julia models/spin1_largeD/test/runtests.jl
#
# It re-verifies, for every results/certificate_D*.json:
#   * the exact integer RP identity -h = sum_a G_a ⊗ G_a and the real on-site term;
#   * that the stored enclosed spectra are byte-identical to what the dense runs wrote
#     (the L=12 spectra are rebuilt from their compact form and their SHA-256 compared with
#     the hash recorded in the certificate);
#   * complete multiplicity accounting (3^L levels) and orthogonality bounds < 1;
#   * outward-rounded log Z intervals and eps_s >= 1 - P_s(n0, beta0),
#     eps_t >= 1 - P_tau(2 n0, beta0) for the stored exact rationals;
#   * the recurrence bound from CertifiedGaps.certified_gap_bound, equal to the stored one.
# It then runs the small dense cross-checks (test/small.jl) and exact symmetry checks
# (test/symmetry.jl).
using Test
include("../seed.jl")
include("../spectra.jl")

"""
    recorded_source_matches(file, hash)

True if `file` (relative to the certificate directory) has SHA-256 `hash`, either as it is
now or with the core package name reverted from `CertifiedGaps` to `ThermalMPOGap`. The dense
runs used the core before its rename; the only later change to the producing scripts is that
literal identifier (the old name never occurs in the current files, so the reversal is exact).
"""
function recorded_source_matches(file, hash)
    text = read(joinpath(root, file), String)
    occursin("ThermalMPOGap", text) && return false
    return text_sha256(text) == hash ||
           text_sha256(replace(text, "CertifiedGaps" => "ThermalMPOGap")) == hash
end

# Scripts that produced the stored data. The recorded hashes of test/*.jl refer to the
# historical checkers, which have since been rewritten, so they are not compared.
const PRODUCERS = ("model.jl", "seed.jl", "certify.jl", "certify_parallel.jl")

"Parse the stored enclosure for (D, L), rebuilding compact files, and return it with its hash."
function stored_enclosure(D, L)
    name = "enclosed_D$(numerator(D))_$(denominator(D))_L$L"
    text = enclosure_text(joinpath(root, "results"), name)
    return JSON3.read(text, Dict{String,Any}), "results/$name.json", text_sha256(text)
end

"All of a spectrum's producing-script hashes that should still match."
producer_hashes(record) = [(f, h) for (f, h) in record["source_hashes"] if f in PRODUCERS]

@testset "Stored spin-one large-D certificates" begin
    names = filter(x -> startswith(x, "certificate_D") && endswith(x, ".json"),
        readdir(joinpath(root, "results")))
    @test length(names) == 2
    for name in names
        data = JSON3.read(read(joinpath(root, "results", name), String), Dict{String,Any})
        D = parse(Rational{BigInt}, data["D"])
        beta = parse(Rational{BigInt}, data["beta"])
        es = parse(Rational{BigInt}, data["eps_s"])
        et = parse(Rational{BigInt}, data["eps_t"])
        n = Int(data["n0"])
        @test rp_identity(D)
        @test iseven(n) && data["certified"]
        @test 0 <= es < 1 && 0 <= et < 1 && beta > 0
        for (file, hash) in producer_hashes(data)
            @test recorded_source_matches(file, hash)
        end
        a, apath, ahash = stored_enclosure(D, n)
        b, bpath, bhash = stored_enclosure(D, 2n)
        @test data["enclosure_hashes"][apath] == ahash
        @test data["enclosure_hashes"][bpath] == bhash
        for spec in (a, b)
            L = spec["L"]
            @test spec["complete"] && spec["rigorous"]
            @test parse(Rational{BigInt}, spec["D"]) == D
            @test sum(s["dimension"] * s["multiplicity"] for s in spec["sectors"]) == 3^L
            @test length(spec["energies"]) == length(spec["radii"]) == 3^L
            @test all(s["levels"] == s["dimension"] for s in spec["sectors"])
            @test all(s["orthogonality_upper"] < 1 for s in spec["sectors"])
            for (file, hash) in producer_hashes(spec)
                @test recorded_source_matches(file, hash)
            end
        end
        za = enclosed_logz(a, beta); zb = enclosed_logz(b, beta); zbb = enclosed_logz(b, 2beta)
        @test es >= Rational{BigInt}(sup(1 - exp(zb - 2za)))    # 1 - Z_2n(b)/Z_n(b)^2
        @test et >= Rational{BigInt}(sup(1 - exp(zbb - 2zb)))   # 1 - Z_2n(2b)/Z_2n(b)^2
        result = certified_gap_bound(es, et, beta)
        @test result.bound !== nothing
        @test result.bound == data["bound"] > 0
        println("D=$D n0=$n beta0=$beta: uniform gap >= $(result.bound) for every even L >= $n")
    end
end

include("small.jl")
include("symmetry.jl")
