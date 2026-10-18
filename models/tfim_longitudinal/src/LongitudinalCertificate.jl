"""
    LongitudinalCertificate

Gap certificates for the longitudinal-field TFIM ring
`H_L = -Σ Z_j Z_{j+1} - g Σ X_j - h Σ Z_j` (Pauli matrices, one site per cell) through the
shared core pipeline. Only the model-specific parts live here: the exact `RingModel` in the
physical σᶻ frame with its lattice reflection and RP factor `C = Z`, and the 30-point grid with
its seeds (n₀, β₀). Certification, checking, the Lean seed table and the summary rows are the
core `certify_point`, `check_point`, `lean_points` and `certificate_row`.
"""
module LongitudinalCertificate
using LinearAlgebra, CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.RingModels
using CertifiedGaps.GapCertificates
using CertifiedGaps.ExactArithmetic: parse_rational
using CertifiedGaps.Spins: exact_pauli
BLAS.set_num_threads(1)
# Cheap environment guard: the certificates are produced and replayed with gradual underflow.
get_zero_subnormals() && error("flush-to-zero is enabled; the enclosures assume gradual underflow")

"The model directory (`models/tfim_longitudinal`); results live in `ROOT/results`."
const ROOT = normpath(joinpath(@__DIR__, ".."))
"The Lean statement whose `points` list states the seed defects of every grid point."
const LEAN_FILE = joinpath(ROOT, "..", "..", "lean", "CertifiedGaps", "Models",
    "IsingLongitudinal.lean")
"Transverse fields g of the certified grid."
const GS = [1//4, 1//2, 3//4, 1//1, 3//2, 2//1]
"Longitudinal fields h of the certified grid."
const HS = [1//20, 1//10, 1//4, 1//2, 1//1]
"Seed length (sites) at every grid point."
const N0 = 6
"""
Seed inverse temperature β₀ per grid point `(g, h)`, chosen by the floating-point scan
(`scripts/scan.jl`: best `certified_gap_bound` over n₀ ∈ {2, 4, 6}, β₀ ∈ {1/4, …, 64}) and kept
fixed because the Lean statement assumes these seeds.
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
by lattice momentum. The site reflection x ↦ −x (on-site P = I) is a symmetry because the bond
`Z ⊗ Z` is mirror symmetric; since the model is real, the core then builds every momentum block
as a real symmetric matrix.
"""
function longitudinal(g, h)
    (; x, z) = exact_pauli()
    RingModel(2, -kron(z, z), -g * x - h * z; reflection = [1 0; 0 1],
        name = "Ising ring, transverse field $g, longitudinal field $h (Pauli)")
end

"Exact RP certificate −h_bond = Z ⊗ Z (one real symmetric factor C = Z, weight 1)."
rp(model) = rp_certificate(model; factors = [exact_pauli().z])

"File tag `g<num>_<den>_h<num>_<den>` of a grid point."
tag(g, h) = "g$(numerator(g))_$(denominator(g))_h$(numerator(h))_$(denominator(h))"
"Path of the certificate of a grid point (the leaderboard links these names)."
certificate_path(g, h) = joinpath(ROOT, "results", "certificate_" * tag(g, h) * ".json.gz")
"Paths of the stored per-point certificates (gzip-compressed JSON)."
result_paths() = sort(filter(p -> startswith(basename(p), "certificate_g") &&
    endswith(p, ".json.gz"), readdir(joinpath(ROOT, "results"); join = true)))

"The keywords that define a grid point for both `certify_point` and `check_point`."
point_definition(g, h) = (; n0 = N0, β0 = BETA0[(g, h)], params = Dict("g" => g, "h" => h),
    model_dir = ROOT)

"""
    certify(g, h; progress = false) -> Dict

Certify the grid point `(g, h)` at its seed (N0, β₀) with the core `certify_point` and write
it to `certificate_path(g, h)`.
"""
function certify(g, h; progress = false)
    model = longitudinal(g, h)
    certify_point(certificate_path(g, h), model; point_definition(g, h)..., rp = rp(model),
        progress)
end

"""
    lean_seeds() -> Dict

The seed defects `(; s, t)` of the Lean statement per grid point `(g, h)`, read by the core
`lean_points`. Throws unless the Lean file lists exactly the 30 grid points, each at the β₀ of
`BETA0`.
"""
function lean_seeds()
    seeds = Dict{Tuple{Rational{BigInt},Rational{BigInt}},Any}()
    for p in lean_points(LEAN_FILE)
        g, h, β0 = p.params
        haskey(BETA0, (g, h)) && BETA0[(g, h)] == β0 ||
            error("Lean point ($g, $h, $β0) is not a grid point at its seed")
        seeds[(g, h)] = (; p.s, p.t)
    end
    length(seeds) == length(points()) || error("the Lean file does not list all grid points")
    return seeds
end

"""
    check(g, h; path = certificate_path(g, h), full = false, lean = lean_seeds()) -> Dict

The core `check_point` of the certificate stored at `path` against the definition of the grid
point `(g, h)` (model, reflection, RP factor, seed, params, model-local source hashes), with
stored defects at most the Lean seeds `lean[(g, h)]`. With `full = true` the seed spectra are
re-enclosed and the fresh defects must not exceed the stored εₛ, εₜ. Throws on any mismatch.
"""
function check(g, h; path = certificate_path(g, h), full = false, lean = lean_seeds())
    model = longitudinal(g, h)
    check_point(path, model; point_definition(g, h)..., rp = rp(model), lean = lean[(g, h)],
        spectra = full)
end

"Check that the stored certificates are exactly the 30 grid points `GS × HS`."
function check_grid(paths)
    sort(basename.(paths)) == sort([basename(certificate_path(g, h)) for (g, h) in points()]) ||
        error("the stored certificates are not exactly the 30 grid points")
    true
end

"""
Exact seed rows from the certificates: the core `certificate_row` with the defects under the
names `epsilon_s`, `epsilon_t` that `lean/scripts/check_seeds.jl` reads, and the file name.
"""
function seed_rows()
    rows = map(result_paths()) do path
        row = certificate_row(read_result(path))
        row["epsilon_s"], row["epsilon_t"] = pop!(row, "eps_s"), pop!(row, "eps_t")
        row["file"] = basename(path)
        row
    end
    sort!(rows; by = r -> (parse_rational(r["g"]), parse_rational(r["h"])))
end

"""
    write_seed_table()

Write `results/seeds.json`, the exact seed data of all certificates in one small file, so that
downstream users (e.g. the Lean seed cross-check) need not decompress the certificates.
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

"""
    check_seed_table()

Require `results/seeds.json` to equal the rows rebuilt from the stored certificates.
"""
function check_seed_table()
    table = read_result(joinpath(ROOT, "results", "seeds.json"))
    table["points"] == parse_exact_json(Results.JSON3.write(seed_rows())) ||
        error("results/seeds.json does not match the certificates")
    true
end

end
