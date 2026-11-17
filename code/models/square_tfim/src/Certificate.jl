"Generation and independent replay of square-torus Ising seed certificates."
module Certificate

using LinearAlgebra, IntervalArithmetic
using CertifiedGaps
using CertifiedGaps.ExactArithmetic: RationalBig, parse_rational, round_up
using CertifiedGaps.Results
using CertifiedGaps.RingModels: encode_model, encode_rp
using CertifiedGaps.VerifiedSpectra

include("SquareIsing.jl")
include("CubeRecurrence.jl")
using .SquareIsing, .CubeRecurrence

export POINTS, build_spectra, certify, check, seed_bounds, exact_strings, certificate_path,
    uniform_bound

const TORI = ((2, 2), (4, 2), (4, 4))
const POINTS = ((g = 18, beta = big(2)//big(9)), (g = 20, beta = big(1)//big(5)),
    (g = 24, beta = big(1)//big(6)))
const FORMAT = "CertifiedGaps square TFIM cube certificate v2"
const SEED_GRID = big(10)^8
const REFINED_BETA = Dict(18 => big(53)//240, 20 => big(157)//800, 24 => big(79)//480)
const REFINED_STEPS = 12
const REFINED_GRID = big(10)^2000

"Encode exact integers and rationals as strings throughout a certificate record."
exact_strings(x::Bool) = x
exact_strings(x::Integer) = string(x)
exact_strings(x::Rational) = string(x)
exact_strings(x::AbstractDict) = Dict(string(k) => exact_strings(v) for (k, v) in x)
exact_strings(x::Union{AbstractArray,Tuple}) = [exact_strings(v) for v in x]
exact_strings(x) = x

"Path of the compressed certificate at field g in a model directory."
certificate_path(root, g) = joinpath(root, "results", "certificate_g$g.json.gz")
manifest_path(root, g) = joinpath(root, "results", "spectra_g$g.json.gz")
token(n, m) = "$(n)x$m"
endpoints(x) = string.([RationalBig(inf(x)), RationalBig(sup(x))])

function model_records(g)
    return Dict(string(m) => exact_strings(Dict("model" => encode_model(row_model(m, g)),
        "rp" => encode_rp(row_rp(m, g)))) for m in (2, 4))
end

"Build and store all three full verified spectra, using the shared spectral implementation."
function build_spectra(root, g)
    records = Dict{String,Any}()
    models = model_records(g)
    for (n, m) in TORI
        seconds = @elapsed spectrum = verified_ring_spectrum(row_model(m, g), n; progress = true)
        sum(length(s.centers) * s.multiplicity for s in spectrum) == 2^(n*m) ||
            error("spectrum does not account for the whole Hilbert space")
        name = "spectrum_g$(g)_$(token(n, m)).json.gz"
        path = joinpath(root, "results", name)
        write_spectrum(path, spectrum)
        records[token(n, m)] = exact_strings(Dict("file" => name,
            "sha256" => sha256_file(path), "summary" => spectrum_summary(spectrum),
            "gap_enclosure" => string.(spectral_gap_enclosure(spectrum)), "seconds" => seconds))
        println("Stored $(token(n, m)): $seconds seconds")
        flush(stdout)
    end
    record = Dict("g" => string(g), "models" => models, "spectra" => records,
        "model_source_sha256" => sha256_file(joinpath(root, "src", "SquareIsing.jl")))
    write_result(manifest_path(root, g), record)
    return record
end

function load_spectra(root, record, g)
    record["g"] == string(g) || error("wrong field in spectrum manifest")
    record["models"] == model_records(g) || error("wrong row model or RP factors")
    record["model_source_sha256"] == sha256_file(joinpath(root, "src", "SquareIsing.jl")) ||
        error("spectra were generated with a different model source")
    spectra = Dict{String,Vector{BlockSpectrum}}()
    for (n, m) in TORI
        r = record["spectra"][token(n, m)]
        name = "spectrum_g$(g)_$(token(n, m)).json.gz"
        r["file"] == name || error("unexpected spectrum path")
        path = joinpath(root, "results", name)
        sha256_file(path) == r["sha256"] || error("spectrum hash mismatch: $path")
        s = read_spectrum(path)
        all(b -> b.multiplicity > 0 && b.radius >= 0 && issorted(b.centers) &&
            all(isfinite, b.centers), s) || error("invalid stored spectral block")
        sum(length(b.centers) * b.multiplicity for b in s) == 2^(n*m) ||
            error("incomplete spectrum")
        exact_strings(spectrum_summary(s)) == r["summary"] || error("spectrum summary mismatch")
        collect(string.(spectral_gap_enclosure(s))) == r["gap_enclosure"] ||
            error("finite gap mismatch")
        spectra[token(n, m)] = s
    end
    return spectra
end

"Enclose the three directional seed defects from the four required log partition sums."
function seed_bounds(spectra, beta)
    a = log_partition_enclosure(spectra["2x2"], beta)
    b = log_partition_enclosure(spectra["4x2"], beta)
    c = log_partition_enclosure(spectra["4x4"], beta)
    d = log_partition_enclosure(spectra["4x4"], 2beta)
    defects = [1-exp(b-2a), 1-exp(c-2b), 1-exp(d-2c)]
    # One extra grid unit makes full replay insensitive to last-place eigensolver variation.
    eps = [max(round_up(RationalBig(sup(x)), SEED_GRID), 0) + 1//SEED_GRID for x in defects]
    return Dict("beta" => string(beta), "n0" => "2", "m0" => "2",
        "log_partition" => Dict("2x2_beta" => endpoints(a), "4x2_beta" => endpoints(b),
            "4x4_beta" => endpoints(c), "4x4_2beta" => endpoints(d)),
        "defect_enclosures" => endpoints.(defects), "defects" => string.(eps))
end

"Exact finite-volume ground-energy, energy-density, and gap enclosures."
function finite_volume_summary(spectra)
    return Dict(token(n, m) => let spectrum = spectra[token(n, m)]
        ground = ground_energy_enclosure(spectrum)
        Dict("sites" => string(n*m),
            "ground_energy" => string.(collect(ground)),
            "ground_energy_per_site" => string.(collect(ground ./ (n*m))),
            "gap" => string.(collect(spectral_gap_enclosure(spectrum))))
    end for (n, m) in TORI)
end

"Repeat the asymmetric recurrence before using the common quadratic tail."
function refined_chain(seed)
    return cube_chain(parse_rational.(seed["defects"]), parse_rational(seed["beta"]);
        min_steps = REFINED_STEPS, grid = REFINED_GRID)
end

"The larger of the two independently checked uniform gap bounds."
uniform_bound(cert) = max(parse_rational(cert["chain"]["gap_lower_bound"]),
    parse_rational(cert["refinement"]["chain"]["gap_lower_bound"]))

"Certify one field and inverse temperature, reusing only model-checked verified spectra."
function certify(root, g, beta; reuse = false)
    record = reuse ? read_result(manifest_path(root, g)) : build_spectra(root, g)
    spectra = load_spectra(root, record, g)
    seed = seed_bounds(spectra, beta)
    chain = cube_chain(parse_rational.(seed["defects"]), beta)
    refined_seed = seed_bounds(spectra, REFINED_BETA[g])
    certificate = Dict("format" => FORMAT, "J" => "1", "g" => string(g),
        "scope" => "unique ground state; every even square torus of side L >= 2",
        "seed" => seed, "chain" => chain, "spectral_record" => record,
        "refinement" => Dict("seed" => refined_seed, "chain" => refined_chain(refined_seed),
            "preliminary_steps" => string(REFINED_STEPS)),
        "finite_volume" => finite_volume_summary(spectra),
        "metadata" => exact_strings(model_metadata(root, Dict("g" => string(g));
            dirs = ("src", "scripts"),
            packages = (CertifiedGaps, IntervalArithmetic), certified = true)))
    write_result(certificate_path(root, g), certificate)
    return certificate
end

"Recompute the complete certificate; with full=true, also rebuild all verified spectra."
function check(root, g; full = false)
    cert = read_result(certificate_path(root, g))
    cert["format"] == FORMAT || error("unsupported certificate format")
    cert["J"] == "1" && cert["g"] == string(g) || error("wrong Hamiltonian parameters")
    cert["scope"] == "unique ground state; every even square torus of side L >= 2" ||
        error("wrong claimed scope")
    verify_model_metadata(cert["metadata"], root; dirs = ("src", "scripts"))
    spectra = load_spectra(root, cert["spectral_record"], g)
    beta = parse_rational(cert["seed"]["beta"])
    seed_bounds(spectra, beta) == cert["seed"] || error("seed replay mismatch")
    eps = parse_rational.(cert["seed"]["defects"])
    cube_chain(eps, beta) == cert["chain"] || error("cube recurrence mismatch")
    refinement = cert["refinement"]
    refinement["preliminary_steps"] == string(REFINED_STEPS) || error("wrong refinement depth")
    seed_bounds(spectra, REFINED_BETA[g]) == refinement["seed"] ||
        error("refined seed replay mismatch")
    refined_chain(refinement["seed"]) == refinement["chain"] ||
        error("refined cube recurrence mismatch")
    finite_volume_summary(spectra) == cert["finite_volume"] ||
        error("finite-volume enclosure mismatch")
    if full
        fresh = Dict(token(n, m) => verified_ring_spectrum(row_model(m, g), n;
            progress = true) for (n, m) in TORI)
        for (n, m) in TORI
            sum(length(s.centers) * s.multiplicity for s in fresh[token(n, m)]) == 2^(n*m) ||
                error("full replay lost spectral multiplicity")
        end
        for stored in (cert["seed"], refinement["seed"])
            replay = seed_bounds(fresh, parse_rational(stored["beta"]))
            limits = parse_rational.(stored["defects"])
            all(parse_rational(replay["defect_enclosures"][i][2]) <= limits[i] for i in 1:3) ||
                error("full spectral replay exceeds stored seed bounds")
        end
    end
    return cert
end

end
