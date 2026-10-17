"""
    LongitudinalCertificate

Gap certificates for the longitudinal-field TFIM ring
`H_L = -Σ Z_j Z_{j+1} - g Σ X_j - h Σ Z_j` (Pauli matrices, one site per cell) through the
shared core pipeline. Only the model-specific parts live here: the exact `RingModel` in the
physical σᶻ frame with its RP factor `C = Z`, the 30-point grid with its seeds (n₀, β₀), and
the loop over `GapCertificates.certify_gap` / `check_certificate`. Block construction,
eigenvalue enclosure, log Z, seed defects and the rational chain are all core code.
"""
module LongitudinalCertificate
using LinearAlgebra, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.RingModels
using CertifiedGaps.GapCertificates: certify_gap, check_certificate
using CertifiedGaps.ExactArithmetic: parse_rational
using CertifiedGaps.Spins: exact_pauli
BLAS.set_num_threads(1)
# Cheap environment guard: the certificates are produced and replayed with gradual underflow.
get_zero_subnormals() && error("flush-to-zero is enabled; the enclosures assume gradual underflow")

"The model directory (`models/tfim_longitudinal`); results live in `ROOT/results`."
const ROOT = normpath(joinpath(@__DIR__, ".."))
"Transverse fields g of the certified grid."
const GS = [1//4, 1//2, 3//4, 1//1, 3//2, 2//1]
"Longitudinal fields h of the certified grid."
const HS = [1//20, 1//10, 1//4, 1//2, 1//1]
"Seed length (sites) at every grid point."
const N0 = 6
"""
Seed inverse temperature β₀ per grid point `(g, h)`, chosen by the floating-point scan
(`scripts/scan.jl`: best `certified_gap_bound` over n₀ ∈ {2, 4, 6}, β₀ ∈ {1/4, …, 64}).
"""
const BETA0 = Dict{Tuple{Rational{Int},Rational{Int}},Rational{Int}}(
    (1//4, 1//20) => 97//4, (1//4, 1//10) => 55//4, (1//4, 1//4) => 13//2,
    (1//4, 1//2) => 11//2, (1//4, 1//1) => 21//4,
    (1//2, 1//20) => 13//1, (1//2, 1//10) => 8//1, (1//2, 1//4) => 9//2,
    (1//2, 1//2) => 4//1, (1//2, 1//1) => 4//1,
    (3//4, 1//20) => 8//1, (3//4, 1//10) => 21//4, (3//4, 1//4) => 7//2,
    (3//4, 1//2) => 7//2, (3//4, 1//1) => 7//2,
    (1//1, 1//20) => 13//2, (1//1, 1//10) => 17//4, (1//1, 1//4) => 3//1,
    (1//1, 1//2) => 3//1, (1//1, 1//1) => 3//1,
    (3//2, 1//20) => 7//2, (3//2, 1//10) => 13//4, (3//2, 1//4) => 11//4,
    (3//2, 1//2) => 5//2, (3//2, 1//1) => 5//2,
    (2//1, 1//20) => 9//4, (2//1, 1//10) => 9//4, (2//1, 1//4) => 9//4,
    (2//1, 1//2) => 9//4, (2//1, 1//1) => 2//1)

"Grid points in g-major order (the order of `certify.jl` point indices)."
points() = [(g, h) for g in GS for h in HS]

"""
    longitudinal(g, h) -> RingModel

The ring in the physical σᶻ frame (the frame unitary is the identity): bond `-Z ⊗ Z`,
on-site term `-g X - h Z`. No charge is conserved for g, h ≠ 0, so the core splits the ring
by lattice momentum only.
"""
function longitudinal(g, h)
    (; x, z) = exact_pauli()
    RingModel(2, -kron(z, z), -g * x - h * z;
        name = "Ising ring, transverse field $g, longitudinal field $h (Pauli)")
end

"Exact RP certificate −h_bond = Z ⊗ Z (one real symmetric factor C = Z, weight 1)."
rp(model) = rp_certificate(model; factors = [exact_pauli().z])

"File tag `g<num>_<den>_h<num>_<den>` of a grid point."
tag(g, h) = "g$(numerator(g))_$(denominator(g))_h$(numerator(h))_$(denominator(h))"
"Path of the core-format certificate of a grid point."
certificate_path(g, h) = joinpath(ROOT, "results", "certificate_" * tag(g, h) * ".json.gz")
"Paths of the per-point certificates (gzip-compressed JSON)."
result_paths() = sort(filter(p -> startswith(basename(p), "certificate_g") &&
    endswith(p, ".json.gz"), readdir(joinpath(ROOT, "results"); join = true)))
"Format an exact rational as a `\"num//den\"` string."
str(x::Rational) = string(numerator(x), "//", denominator(x))

"""
Model-local sources that produce and check the certificates, hashed into each result's
`metadata` by the core `model_metadata` and enforced by `verify_model_metadata`. The
environment is the core's (`julia/Project.toml`, `julia/Manifest.toml`), recorded as
provenance under `core_sha256`.
"""
const PROVENANCE = (dirs = ("src",), files = ("certify.jl", "check.jl"))

"Core `model_metadata` (model-local `source_sha256`, core `core_sha256`) for a result."
metadata(config; certified = false) =
    model_metadata(ROOT, config; PROVENANCE..., packages = (CertifiedGaps,), certified)

"""
    certify_point(g, h; progress = false) -> Dict

`certify_gap` at the grid point's seed (N0, β₀), with the fields and run metadata attached.
"""
function certify_point(g, h; progress = false)
    model = longitudinal(g, h)
    cert = certify_gap(model; n0 = N0, β0 = BETA0[(g, h)], rp = rp(model), progress)
    cert["g"], cert["h"] = str(g), str(h)
    cert["metadata"] = metadata(Dict("g" => str(g), "h" => str(h)); certified = true)
    cert
end

"""
    check_point(path; full = false) -> Dict

`check_record` of the certificate stored at `path`, named by its file.
"""
check_point(path; full = false) = check_record(read_result(path), basename(path); full)

"""
    check_record(cert, name; full = false) -> Dict

Verify one certificate (as read by `read_result`, stored under the file name `name`) and return
the recomputed values; throw on any mismatch. Enforces the recorded model-local source hashes
(and that their key set is exactly the current model-local set), requires the stored model, RP
data and seed to be exactly those of the grid point, then runs the core `check_certificate`:
with `full = false` the exact RP identity and the exact rational chain from the stored seed;
with `full = true` also freshly verified spectra of the rings of 6 and 12 sites, whose defect
enclosures must not exceed the stored εₛ, εₜ.
"""
function check_record(cert, name; full = false)
    verify_model_metadata(cert["metadata"], ROOT; PROVENANCE...)
    g, h = parse_rational(cert["g"]), parse_rational(cert["h"])
    haskey(BETA0, (g, h)) || error("$name: ($g, $h) is not a grid point")
    name == basename(certificate_path(g, h)) || error("$name: misnamed")
    model, expected = decode_model(cert["model"]), longitudinal(g, h)
    model.bond == expected.bond && model.onsite == expected.onsite ||
        error("$name: stored model differs from the longitudinal TFIM at ($g, $h)")
    decode_rp(cert["reflection_positivity"]).factors == rp(expected).factors ||
        error("$name: stored RP factors differ from C = Z")
    seed = cert["seed"]
    parse(Int, seed["n0"]) == N0 && parse_rational(seed["beta0"]) == BETA0[(g, h)] ||
        error("$name: stored seed differs from the grid table")
    check_certificate(cert; spectra = full)
end

"Check that the stored certificates are exactly the 30 grid points `GS × HS`."
function check_grid(paths)
    sort(basename.(paths)) == sort([basename(certificate_path(g, h)) for (g, h) in points()]) ||
        error("the stored certificates are not exactly the 30 grid points")
    true
end

"Exact seed rows (g, h, n₀, β₀, εₛ, εₜ, chain and interval bounds) from the certificates."
function seed_rows()
    rows = Dict{String,Any}[]
    for path in result_paths()
        c = read_result(path); s = c["seed"]
        push!(rows, Dict{String,Any}("g" => c["g"], "h" => c["h"], "n0" => s["n0"],
            "beta0" => s["beta0"], "epsilon_s" => s["eps_s"], "epsilon_t" => s["eps_t"],
            "gap_lower_bound" => c["gap_lower_bound"],
            "interval_bound_decimal" => c["interval_bound_decimal"],
            "file" => basename(path)))
    end
    sort!(rows; by = r -> (parse_rational(r["g"]), parse_rational(r["h"])))
end

"""
    write_seed_table()

Write `results/seeds.json`, the exact seed data of all certificates in one small file, so that
downstream users (e.g. the Lean statement) need not decompress the certificates.
"""
function write_seed_table()
    write_result(joinpath(ROOT, "results", "seeds.json"), Dict(
        "hamiltonian" => "H_L = -sum_j Z_j Z_{j+1} - g sum_j X_j - h sum_j Z_j (Pauli, periodic)",
        "rp" => "-h_bond = Z (x) Z with C = Z real; onsite o = -g X - h Z real symmetric",
        "epsilon_s" => "upper bound on 1 - Z_{2n0}(beta0)/Z_{n0}(beta0)^2",
        "epsilon_t" => "upper bound on 1 - Z_{2n0}(2 beta0)/Z_{2n0}(beta0)^2",
        "gap_lower_bound" => "exact rational chain: gap > this for every even L >= n0 sites",
        "points" => seed_rows()))
end

end
