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
const LEGACY_POINTS = ((g = 18, beta = big(2)//big(9)), (g = 20, beta = big(1)//big(5)),
    (g = 24, beta = big(1)//big(6)))
const FORMAT = "CertifiedGaps square TFIM cube certificate v3"
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

const POINTS = ((g = 15, beta = big(47)//200), (g = 16, beta = big(141)//640),
    (g = 18, beta = big(143)//720), (g = 20, beta = big(29)//160),
    (g = 24, beta = big(5)//32))
const DIRECT_GRID = big(10)^12
const DIRECT_CHAIN_GRID = big(10)^8000
const PARAMETER_RADIUS = big(1)//1000
const SCOPE = "unique ground state; every even square torus of side L >= 4"

"Read quoted exact numbers without regex matching the long rational strings in the tail."
read_certificate(root, g) = read_result(certificate_path(root, g), Dict{String,Any})

# Reject unquoted scientific numbers instead of allowing JSON's Float64 conversion to matter.
function require_exact_strings(x)
    if x isa AbstractDict
        foreach(require_exact_strings, values(x))
    elseif x isa Union{AbstractArray,Tuple}
        foreach(require_exact_strings, x)
    elseif x isa Number && !(x isa Bool)
        error("scientific numbers must be quoted exact strings")
    end
    return nothing
end

"Enclose defects of Py(4,2,2β) and Pτ(4,4,β) directly from complete verified spectra."
function direct_input_bounds(spectra, beta)
    beta = RationalBig(beta)
    beta > 0 || error("inverse temperature must be positive")
    a = log_partition_enclosure(spectra["4x2"], 2beta)
    b = log_partition_enclosure(spectra["4x4"], 2beta)
    c = log_partition_enclosure(spectra["4x4"], beta)
    defects = [1-exp(b-2a), 1-exp(b-2c)]
    eps = [max(round_up(RationalBig(sup(x)), DIRECT_GRID), 0) + 1//DIRECT_GRID
        for x in defects]
    return Dict("beta" => string(beta),
        "log_partition" => Dict("4x2_2beta" => endpoints(a),
            "4x4_2beta" => endpoints(b), "4x4_beta" => endpoints(c)),
        "defect_enclosures" => endpoints.(defects), "defects" => string.(eps))
end

"Apply the direct-temperature prefix and the exact chain at spatial seed side four."
function direct_seed_chain(input_defects, beta)
    beta = RationalBig(beta)
    beta > 0 || error("inverse temperature must be positive")
    defects = direct_temperature_prefix(input_defects...; grid = DIRECT_GRID)
    seed = Dict("n0" => "4", "m0" => "4", "beta" => string(2beta),
        "defects" => string.(collect(defects)), "rounding_grid" => string(DIRECT_GRID))
    chain = cube_chain(defects, 2beta; min_steps = REFINED_STEPS,
        grid = DIRECT_CHAIN_GRID, min_even_side = 4)
    return Dict{String,Any}("seed" => seed, "chain" => chain)
end

"A direct-temperature certificate, with the spectral inputs separated from its exact prefix."
function direct_temperature_record(spectra, beta)
    input = direct_input_bounds(spectra, beta)
    record = direct_seed_chain(parse_rational.(input["defects"]), beta)
    record["input"] = input
    return record
end

"""
    parameter_interval_record(g, beta, input_defects; radius = 1//1000)

Certify all fields in [g-radius,g+radius]. For K=H/g and fixed scaled temperature s=2βg,
Weyl bounds give a loss at most 64sr in each input log purity, where r bounds |1/g'-1/g|.
Use exp(-64sr) >= 1-64sr, then rescale the resulting gap by the smallest field ratio.
"""
function parameter_interval_record(g, beta, input_defects; radius = PARAMETER_RADIUS)
    g, beta, radius = RationalBig(g), RationalBig(beta), RationalBig(radius)
    beta > 0 && 0 < radius < g || error("invalid field interval or temperature")
    lower, upper = g-radius, g+radius
    inverse_radius = radius / (g*lower)
    scaled_beta = 2beta*g
    loss = 64scaled_beta*inverse_radius
    0 < loss < 1 || error("parameter loss leaves the admissible region")
    all(x -> 0 <= x < 1, input_defects) || error("invalid parameter input defects")
    defects = [1-(1-x)*(1-loss) for x in RationalBig.(input_defects)]
    record = direct_seed_chain(defects, beta)
    record["field_interval"] = string.([lower, upper])
    record["inverse_field_radius"] = string(inverse_radius)
    record["scaled_effective_beta"] = string(scaled_beta)
    record["log_purity_loss"] = string(loss)
    record["input_defects"] = string.(defects)
    record["gap_lower_bound"] = string(lower/g *
        parse_rational(record["chain"]["gap_lower_bound"]))
    record["scope"] = SCOPE * "; every field in the closed interval; J = 1"
    return record
end

"Preserve the earlier point certificates without duplicating their spectral evidence."
function legacy_record(spectra, g)
    p = findfirst(p -> p.g == g, LEGACY_POINTS)
    isnothing(p) && return nothing
    beta = LEGACY_POINTS[p].beta
    seed = seed_bounds(spectra, beta)
    refined_seed = seed_bounds(spectra, REFINED_BETA[g])
    return Dict("seed" => seed,
        "chain" => cube_chain(parse_rational.(seed["defects"]), beta),
        "refinement" => Dict("seed" => refined_seed, "chain" => refined_chain(refined_seed),
            "preliminary_steps" => string(REFINED_STEPS)))
end

"The strongest checked point bound, always with the primary certificate's L >= 4 scope."
function uniform_bound(cert)
    direct = parse_rational(cert["direct_temperature"]["chain"]["gap_lower_bound"])
    legacy = cert["legacy"]
    isnothing(legacy) && return direct
    return max(direct, parse_rational(legacy["chain"]["gap_lower_bound"]),
        parse_rational(legacy["refinement"]["chain"]["gap_lower_bound"]))
end

"Build all scientific fields through one shared path for generation and replay."
function certificate_payload(record, spectra, g, beta)
    direct = direct_temperature_record(spectra, beta)
    return Dict("format" => FORMAT, "J" => "1", "g" => string(g), "scope" => SCOPE,
        "direct_temperature" => direct,
        "parameter_interval" => parameter_interval_record(g, beta,
            parse_rational.(direct["input"]["defects"])),
        "legacy" => legacy_record(spectra, g), "spectral_record" => record,
        "finite_volume" => finite_volume_summary(spectra))
end

"Certify one configured field, preserving prior bounds and adding a closed field interval."
function certify(root, g, beta; reuse = false)
    any(p -> p.g == g, POINTS) || error("field is not a configured certificate point")
    record = reuse ? read_result(manifest_path(root, g)) : build_spectra(root, g)
    spectra = load_spectra(root, record, g)
    certificate = certificate_payload(record, spectra, g, RationalBig(beta))
    certificate["metadata"] = exact_strings(model_metadata(root, Dict("g" => string(g));
        dirs = ("src", "scripts"), packages = (CertifiedGaps, IntervalArithmetic),
        certified = true))
    write_result(certificate_path(root, g), certificate)
    return certificate
end

"Replay a supplied record; full=true independently regenerates all verified spectral blocks."
function check_record(root, g, cert; full = false)
    any(p -> p.g == g, POINTS) || error("field is not a configured certificate point")
    cert["format"] == FORMAT || error("unsupported certificate format")
    cert["J"] == "1" && cert["g"] == string(g) || error("wrong Hamiltonian parameters")
    cert["scope"] == SCOPE || error("wrong claimed scope")
    for key in ("direct_temperature", "parameter_interval", "legacy", "finite_volume")
        require_exact_strings(cert[key])
    end
    verify_model_metadata(cert["metadata"], root; dirs = ("src", "scripts"))
    spectra = load_spectra(root, cert["spectral_record"], g)
    beta = parse_rational(cert["direct_temperature"]["input"]["beta"])
    expected = certificate_payload(cert["spectral_record"], spectra, g, beta)
    Set(keys(cert)) == union(Set(keys(expected)), Set(["metadata"])) ||
        error("unexpected certificate fields")
    for (key, value) in expected
        cert[key] == value || error("certificate replay mismatch: $key")
    end
    if full
        fresh = Dict(token(n, m) => verified_ring_spectrum(row_model(m, g), n;
            progress = true) for (n, m) in TORI)
        for (n, m) in TORI
            sum(length(s.centers) * s.multiplicity for s in fresh[token(n, m)]) == 2^(n*m) ||
                error("full replay lost spectral multiplicity")
        end
        replay = direct_input_bounds(fresh, beta)
        limits = parse_rational.(cert["direct_temperature"]["input"]["defects"])
        all(parse_rational(replay["defect_enclosures"][i][2]) <= limits[i] for i in 1:2) ||
            error("full spectral replay exceeds direct-temperature input bounds")
        if !isnothing(cert["legacy"])
            legacy = cert["legacy"]
            for stored in (legacy["seed"], legacy["refinement"]["seed"])
                replay = seed_bounds(fresh, parse_rational(stored["beta"]))
                limits = parse_rational.(stored["defects"])
                all(parse_rational(replay["defect_enclosures"][i][2]) <= limits[i] for i in 1:3) ||
                    error("full spectral replay exceeds legacy seed bounds")
            end
        end
    end
    return cert
end

"Recompute a stored certificate, optionally including a complete fresh spectral replay."
check(root, g; full = false) = check_record(root, g, read_certificate(root, g); full)

end
