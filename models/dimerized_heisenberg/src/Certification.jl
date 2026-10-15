# Seed certificates, all-beta exclusions and their checks for the dimerised Heisenberg ring.
#
# `generate(delta)` writes results/seed_delta_<p>_<q>.json.gz: exact RP factors, complete
# verified 4/8/16-spin spectra, the enclosed seeds at n0 = 2, 4 cells (or a proof that every
# beta > 0 fails), the survey seed, the 16-spin gap interval and the dimer perturbation bound.
# `check_record(record)` recomputes everything from fresh spectra (check.jl --full);
# `check_record(record; spectra=stored_spectra(record))` uses the stored centres after
# re-deriving every block's radius from its stored witness (check.jl, fast).
module Certification

using CertifiedGaps, CertifiedGaps.Results, IntervalArithmetic, LinearAlgebra
using CertifiedGaps.ExactArithmetic: round_up, round_down
using CertifiedGaps.RingModels: RingModel, rp_certificate, verify_rp, encode_model,
    encode_rp, decode_model, decode_rp
include("verified_spectra.jl")
include("rp.jl")
using .VerifiedSpectra, .DimerizedRP

export generate, check_record, load_record, seed_intervals, dimer_bounds, provenance,
    core_provenance, exact, endpoints, DELTAS, negative_cover, replay_spectra, ring_levels, stored_spectra,
    seed_path

const Q = Rational{BigInt}
const MODEL_DIR = normpath(joinpath(@__DIR__, ".."))
const ROOT = normpath(joinpath(MODEL_DIR, "../.."))
const DELTAS = Q[1//2, 3//10, 1//4, 1//5, 3//20, 1//10]
const FORMAT = "dimerized-heisenberg-v2"
exact(s) = (p = split(s, "//"); parse(BigInt, p[1]) // parse(BigInt, p[2]))
endpoints(x) = string.([Q(inf(x)), Q(sup(x))])
qinterval(x) = CertifiedGaps.to_interval(x)

"""
    provenance()

SHA-256 of the model-local sources that produce or check a seed certificate (by repository
path). `check_record` requires exactly these keys and every digest to match. The shared core
is recorded separately by `core_provenance` and is not enforced.
"""
function provenance()
    m = "models/dimerized_heisenberg"
    source_hashes(ROOT; dirs=(), files=(
        "$m/src/Certification.jl", "$m/src/rp.jl", "$m/src/sectors.jl",
        "$m/src/verified_spectra.jl", "$m/src/provenance.jl", "$m/certify.jl",
        "$m/check.jl", "$m/Project.toml", "$m/Manifest.toml"))
end

"""
    core_provenance()

SHA-256 of the shared core (`julia/src`, `julia/Project.toml`, `julia/Manifest.toml`), stored
as `metadata["core_sha256"]`. Provenance only: no checker compares it with the checkout, so a
core edit does not invalidate a certificate (the checker reruns the core arithmetic instead).
"""
core_provenance() = source_hashes(ROOT; dirs=("julia/src",),
    files=("julia/Project.toml", "julia/Manifest.toml"))

"Path of the seed certificate of `delta`."
seed_path(delta) = joinpath(MODEL_DIR, "results",
    "seed_delta_$(numerator(delta))_$(denominator(delta)).json.gz")

"Load exact JSON, including gzip artifacts, without routing integer tokens through Float64."
load_record(path) = parse_exact_json(endswith(path, ".gz") ?
    read(`gzip -dc $path`, String) : read(path, String))
asint(x) = x isa AbstractString ? parse(Int, x) : Int(x)

"Recover all eigenvalue intervals with every sector multiplicity."
function ring_levels(ring)
    energies = Float64[]
    radii = Q[]
    for block in ring["blocks"]
        vals = exact.(block["eigenvalues"])
        length(vals) == asint(block["dimension"]) || error("block dimension mismatch")
        all(x -> Q(Float64(x)) == x, vals) || error("nonbinary center")
        for _ in 1:asint(block["multiplicity"])
            append!(energies, Float64.(vals))
            append!(radii, fill(exact(block["eigenvalue_radius"]), length(vals)))
        end
    end
    length(energies) == 2^asint(ring["spins"]) || error("incomplete spectrum")
    return energies, radii
end

"""
    stored_spectra(record)

Eigenvalue centres and radii of the stored rings, for the fast check. Every block's residual,
orthogonality and radius witness is re-derived from its stored component maxima with
`scalar_witness` and must equal the stored values; the sector dimensions must exhaust every
charge; a block reused by momentum conjugation must repeat its partner exactly. This trusts the
stored maxima of the eigenvector products; `replay_spectra` (check.jl --full) redoes them.
"""
function stored_spectra(record)
    spectra = Dict{Int,Any}()
    parity(b) = b["spin_flip_parity"] === nothing ? 0 : asint(b["spin_flip_parity"])
    for ring in record["rings"]
        L = asint(ring["spins"])
        bylabel = Dict{Tuple{Int,Int,Int},Any}()
        for b in ring["blocks"]
            w = scalar_witness(asint(b["dimension"]), L, exact(b["h_component_max"]),
                exact(b["q_component_max"]), exact(b["lambda_max_abs"]),
                exact(b["raw_residual_component_max"]),
                exact(b["raw_orthogonality_component_max"]), exact(b["matrix_error_upper"]))
            w.steps * VerifiedSpectra.UNIT < 1 || error("rounding model out of range")
            w.coarse_ebound == exact(b["coarse_orthogonality_norm_upper"]) < 1 ||
                error("coarse orthogonality witness mismatch")
            w.ebound == exact(b["orthogonality_norm_upper"]) < 1 ||
                error("orthogonality witness mismatch")
            w.rbound == exact(b["residual_norm_upper"]) || error("residual witness mismatch")
            w.radius == exact(b["eigenvalue_radius"]) > 0 || error("radius witness mismatch")
            issorted(exact.(b["eigenvalues"])) || error("unsorted block spectrum")
            label = (asint(b["down_spins"]), asint(b["cell_momentum"]), parity(b))
            haskey(bylabel, label) && error("duplicate sector $label")
            bylabel[label] = b
        end
        n = L ÷ 2
        for ((q, k, p), b) in bylabel
            haskey(b, "conjugate_momentum") || continue
            m = asint(b["conjugate_momentum"])
            m == n - k && k > n ÷ 2 || error("invalid conjugate momentum")
            partner = bylabel[(q, m, p)]
            haskey(partner, "conjugate_momentum") && error("conjugate of a conjugate")
            for key in ("dimension", "multiplicity", "eigenvalues", "eigenvalue_radius")
                b[key] == partner[key] || error("conjugate block differs in $key")
            end
        end
        for q in 0:L÷2
            sum(asint(b["dimension"]) for b in ring["blocks"]
                if asint(b["down_spins"]) == q; init=0) == binomial(L, q) ||
                error("charge $q incomplete")
        end
        spectra[L] = ring_levels(ring)
    end
    return spectra
end

"All three log-partition intervals and both purity defects for an even seed."
function seed_intervals(spectra, n, beta)
    iseven(n) || error("odd seeds are outside the stated theorem")
    small = enclosed_logZ(spectra[2n]..., beta)
    large = enclosed_logZ(spectra[4n]..., beta)
    double = enclosed_logZ(spectra[4n]..., 2beta)
    return (; s=1-exp(large-2small), t=1-exp(double-2large),
        logs=Dict("n_beta"=>endpoints(small), "2n_beta"=>endpoints(large),
            "2n_2beta"=>endpoints(double)))
end

"Exact scalable norm-perturbation route about disconnected strong dimers."
function dimer_bounds(delta::Q, n::Int, beta::Q)
    # ||V_n|| <= 3 n (1-delta)/4; each log Z error is <= beta ||V_n||.
    v = 3n*(1-delta)/4
    x = exp(-qinterval(beta*(1+delta)))
    single = (1+3x^2)/(1+3x)^2
    ps_lower = exp(-qinterval(4beta*v))
    pt_lower = single^(2n)*exp(-qinterval(8beta*v))
    return (; s=1-ps_lower, t=1-pt_lower, norm_small=v)
end

"Reduced partition sum with an arbitrary exact energy shift, for a beta interval."
function reduced_logZ(E, r, beta)
    shift = minimum(E)
    total = interval(0.0)
    cache = Dict(radius=>sup(qinterval(radius)) for radius in unique(r))
    for (e, radius) in zip(E, r)
        rr = cache[radius]
        energy = interval(e)-interval(shift)+interval(-rr, rr)
        total += exp(-beta*energy)
    end
    return log(total), Q(shift)
end

"Lower bounds on both defects throughout a closed beta interval."
function lower_defects(small, large, a::Q, b::Q)
    beta = interval(inf(qinterval(a)), sup(qinterval(b)))
    zn, en = reduced_logZ(small..., beta)
    z2n, e2n = reduced_logZ(large..., beta)
    z2n2, _ = reduced_logZ(large..., 2beta)
    ps = exp(-beta*qinterval(e2n-2en)+z2n-2zn)
    pt = exp(z2n2-2z2n)
    return max(0//big(1), Q(inf(1-ps))), max(0//big(1), Q(inf(1-pt)))
end

"Prove that even optimistic defects fail the positive-rate doubling test."
function rejects(s::Q, t::Q; steps=16)
    # Round downward after every step; monotonicity preserves the optimistic lower bounds.
    grid = big(10)^15
    for j in 0:steps
        z = 1-(1-s)^2*(1-t)
        z >= 1//2 && return j
        s = round_down(z^2/(2*(1-z)^2), grid)
        zt = 1-(1-t)^2*(1-s)
        zt >= 1//2 && return j+1  # Next z(s',t') >= zt^2/(2(1-zt)^2) >= 1/2.
        t = round_down(zt^2/(2*(1-zt)^2), grid)
    end
    return nothing
end

"Cover all beta>0 by rigorously rejected intervals, or return an unresolved interval."
function negative_cover(spectra, n; max_depth=16)
    small, large = spectra[2n], spectra[4n]
    # Thermal purity is nondecreasing in beta. If P_tau(2n,1)<=1/2, all beta<=1 fail.
    first = seed_intervals(spectra, n, 1//big(1))
    inf(first.t) >= 1//2 || return Dict("status"=>"unresolved_low_beta")
    # E0(2n)-2E0(n)>0 and the trivial dimension bound imply
    # Ps <= dim(2n) exp[-beta (E0(2n)-2E0(n))].
    en, rn = small
    e2, r2 = large
    elo2 = minimum(Q(e)-r for (e,r) in zip(e2,r2))
    ehin = minimum(Q(e)+r for (e,r) in zip(en,rn))
    diff = elo2-2ehin
    diff > 0 || return Dict("status"=>"unresolved_high_beta")
    tail = big(32)//1
    while sup(exp(qinterval(4n)*log(interval(2))-qinterval(tail*diff))) > 1//2
        tail *= 2
    end
    cover = Any[]
    unresolved = Any[]
    function visit(a, b, depth)
        s, t = lower_defects(small, large, a, b)
        j = rejects(s, t)
        if j !== nothing
            push!(cover, Dict("a"=>string(a), "b"=>string(b), "s_lower"=>string(s),
                "t_lower"=>string(t), "rejected_by_step"=>string(j)))
        elseif depth == max_depth
            push!(unresolved, [string(a), string(b)])
        else
            mid = (a+b)/2
            visit(a, mid, depth+1)
            visit(mid, b, depth+1)
        end
    end
    visit(big(1)//1, tail, 0)
    return Dict("status"=>isempty(unresolved) ? "all_beta_rejected" : "unresolved",
        "n0_cells"=>string(n), "low_beta_end"=>"1//1", "tail_beta"=>string(tail),
        "ground_energy_defect_lower"=>string(diff), "cover"=>cover,
        "unresolved"=>unresolved,
        "scope"=>"No positive rate from the stated doubling recurrence, not a gaplessness claim")
end

"Float64 grid used only to choose candidates; every positive claim is then enclosed."
function scan_candidates(spectra)
    rows = Any[]
    specs = Dict(L=>VerifiedSpectra.PartialSpectrum(sort(E), 2^L, true)
        for (L,(E,r)) in spectra)
    for n in (2, 4), bi in 1:128
        beta = bi//big(4)
        small, large = specs[2n], specs[4n]
        zn = VerifiedSpectra.log_partition(small.E, Float64(beta))
        zl = VerifiedSpectra.log_partition(large.E, Float64(beta))
        zd = VerifiedSpectra.log_partition(large.E, 2Float64(beta))
        s, t = -expm1(zl-2zn), -expm1(zd-2zl)
        bound = s < 1 && t < 1 ?
            certified_gap_bound(max(s,0.0), max(t,0.0), beta).bound : nothing
        push!(rows, Dict("n0_cells"=>string(n), "beta0"=>string(beta),
            "epsilon_s_float"=>s, "epsilon_t_float"=>t, "bound_float"=>bound))
    end
    return rows
end

function positive_seed(spectra, n, beta)
    fresh = seed_intervals(spectra, n, beta)
    # A decimal grid with a spare unit makes independent BLAS replays portable.
    es = round_up(Q(sup(fresh.s)), big(10)^6)+1//big(10)^6
    et = round_up(Q(sup(fresh.t)), big(10)^6)+1//big(10)^6
    result = certified_gap_bound(es, et, beta)
    result.bound === nothing && return nothing
    bound = round_down(Q(result.bound), big(10)^6)
    logs = Dict(key=>string.([round_down(exact(val[1]),big(10)^6)-1//big(10)^6,
        round_up(exact(val[2]),big(10)^6)+1//big(10)^6]) for (key,val) in fresh.logs)
    return Dict("n0_cells"=>string(n), "beta0"=>string(beta),
        "eps_s"=>string(es), "eps_t"=>string(et), "gap_lower"=>string(bound),
        "log_partition"=>logs,
        "validity_sites"=>"N divisible by 4, N >= $(2n)",
        "interval_bound_binary"=>string(Q(result.bound)))
end

"Generate current RP, verified spectra, seed and documented negative results for one delta."
function generate(delta::Q; exclude_all_beta=true)
    delta in DELTAS || error("parameter outside requested grid")
    BLAS.set_num_threads(1)
    start = time()
    rp = audit_rp(delta)
    rp.conjugate_identity || error("RP identity failed")
    model = RingModel(4,rp.transformed_bond,rp.transformed_onsite)
    factors = rp_certificate(model; factors=rp.channels,weights=fill(1-delta,3))
    spectra = Dict{Int,Any}()
    rings = Any[]
    for L in (4,8,16)
        ring, E, r = verified_ring(delta, L)
        push!(rings, ring)
        spectra[L] = (E,r)
    end
    scan = scan_candidates(spectra)
    candidates = filter(r -> r["bound_float"] !== nothing, scan)
    seeds = Any[]
    for n in (2,4)
        eligible = filter(r -> r["n0_cells"] == string(n), candidates)
        isempty(eligible) && continue
        best = argmax(r -> r["bound_float"], eligible)
        seed = positive_seed(spectra, n, exact(best["beta0"]))
        seed === nothing || push!(seeds, seed)
    end
    # Re-certify the survey parameters independently of the beta-grid optimization.
    survey_beta = delta == 1//2 ? 25//big(4) : delta == 3//10 ? 13//big(2) : nothing
    survey = survey_beta === nothing ? nothing : positive_seed(spectra, 4, survey_beta)
    exclusions = Any[]
    if exclude_all_beta
        for n in (2,4)
            any(s -> s["n0_cells"] == string(n), seeds) && continue
            push!(exclusions, negative_cover(spectra, n))
        end
    end
    perturbation = dimer_bounds(delta, 4, survey_beta === nothing ? 8//big(1) : survey_beta)
    E, r = spectra[16]
    lows = sort([Q(e)-rr for (e,rr) in zip(E,r)])
    highs = sort([Q(e)+rr for (e,rr) in zip(E,r)])
    gap = string.([round_down(lows[2]-highs[1],big(10)^6)-1//big(10)^6,
        round_up(highs[2]-lows[1],big(10)^6)+1//big(10)^6])
    record = Dict("format"=>FORMAT, "delta"=>string(delta), "rings"=>rings,
        "model"=>encode_model(model), "reflection_positivity"=>encode_rp(factors),
        "rp_identity"=>"-h=(1-delta) sum C tensor conj(C), C=M' B M/2",
        "seeds"=>seeds, "survey_seed"=>survey, "negative_results"=>exclusions,
        "scan_float_not_certificate"=>scan, "reference_gap_16_interval"=>gap,
        "dimer_perturbation"=>Dict("n0"=>"4",
            "beta0"=>string(survey_beta === nothing ? 8//big(1) : survey_beta),
            "norm_small"=>string(perturbation.norm_small),
            "eps_s_enclosure"=>endpoints(perturbation.s),
            "eps_t_enclosure"=>endpoints(perturbation.t)),
        "seconds_total"=>time()-start,
        "metadata"=>merge(run_metadata(Dict("delta"=>string(delta));
                source_sha256=provenance(), packages=(CertifiedGaps,IntervalArithmetic),
                status="recertified_v2", certified=!isempty(seeds)),
            Dict("core_sha256"=>core_provenance())))
    path = chopsuffix(seed_path(delta), ".gz")
    write_result(path, record)
    run(`gzip -n -f $path`)
    println("SAVED $(path).gz seconds=$(record["seconds_total"]) seeds=$seeds")
    return record
end

"Rebuild every Hamiltonian block and residual; saved residual maxima are never trusted."
function replay_spectra(record)
    delta = exact(record["delta"])
    spectra = Dict{Int,Any}()
    for ring in record["rings"]
        L = asint(ring["spins"])
        rebuilt, E, r = verified_ring(delta, L)
        # Check complete labels and dimensions, not just the total dimension.
        labels(b) = (asint(b["down_spins"]), asint(b["cell_momentum"]),
            b["spin_flip_parity"] === nothing ? 0 : asint(b["spin_flip_parity"]))
        old = Dict(labels(b)=>b for b in ring["blocks"])
        length(old) == length(rebuilt["blocks"]) || error("sector list mismatch")
        for b in rebuilt["blocks"]
            stored = old[labels(b)]
            for key in ("dimension", "multiplicity")
                asint(b[key]) == asint(stored[key]) || error("sector $key mismatch")
            end
            vals, oldvals = exact.(b["eigenvalues"]), exact.(stored["eigenvalues"])
            rho = exact(b["eigenvalue_radius"])+exact(stored["eigenvalue_radius"])
            all(abs.(vals-oldvals) .<= rho) || error("stored spectrum disagrees with fresh solve")
        end
        spectra[L] = (E,r)
    end
    return spectra
end

"Check the actual theorem inputs with fresh full-spectrum verification."
function check_record(record; spectra=nothing, check_hashes=true)
    record["format"] == FORMAT || error("wrong certificate format")
    classified = vcat([asint(s["n0_cells"]) for s in record["seeds"]],
        [asint(r["n0_cells"]) for r in record["negative_results"]])
    sort(classified) == [2,4] || error("each eligible seed must be classified exactly once")
    if check_hashes
        Set(keys(record["metadata"]["source_sha256"])) == Set(keys(provenance())) ||
            error("source fingerprint list mismatch")
        for (path, digest) in record["metadata"]["source_sha256"]
            sha256_file(joinpath(ROOT,path)) == digest || error("source hash mismatch: $path")
        end
    end
    rp = audit_rp(exact(record["delta"]))
    model = RingModel(4,rp.transformed_bond,rp.transformed_onsite)
    encode_model(model) == encode_model(decode_model(record["model"])) ||
        error("stored physical model mismatch")
    verify_rp(model,decode_rp(record["reflection_positivity"])) || error("stored RP failed")
    spectra === nothing && (spectra = replay_spectra(record))
    seeds = copy(record["seeds"])
    record["survey_seed"] === nothing || push!(seeds, record["survey_seed"])
    for seed in seeds
        n, beta = asint(seed["n0_cells"]), exact(seed["beta0"])
        fresh = seed_intervals(spectra, n, beta)
        es, et = exact(seed["eps_s"]), exact(seed["eps_t"])
        Q(sup(fresh.s)) <= es || error("spatial seed fails fresh verification")
        Q(sup(fresh.t)) <= et || error("thermal seed fails fresh verification")
        bound = certified_gap_bound(es, et, beta).bound
        bound !== nothing && exact(seed["gap_lower"]) <= Q(bound) || error("gap bound fails")
        Q(bound) == exact(seed["interval_bound_binary"]) || error("recurrence mismatch")
        for (key,values) in fresh.logs
            stored = exact.(seed["log_partition"][key])
            lo,hi = exact.(values)
            stored[1] <= lo <= hi <= stored[2] || error("log-partition enclosure mismatch")
        end
        seed["validity_sites"] == "N divisible by 4, N >= $(2n)" || error("invalid size claim")
    end
    # The dimer perturbation record is a documented (failing) route; recompute it exactly.
    dimer = record["dimer_perturbation"]
    d = dimer_bounds(exact(record["delta"]), asint(dimer["n0"]), exact(dimer["beta0"]))
    dimer["norm_small"] == string(d.norm_small) && dimer["eps_s_enclosure"] == endpoints(d.s) &&
        dimer["eps_t_enclosure"] == endpoints(d.t) || error("dimer perturbation mismatch")
    E,r = spectra[16]
    lows = sort([Q(e)-rr for (e,rr) in zip(E,r)])
    highs = sort([Q(e)+rr for (e,rr) in zip(E,r)])
    lo,hi = exact.(record["reference_gap_16_interval"])
    lo <= lows[2]-highs[1] <= highs[2]-lows[1] <= hi || error("reference gap mismatch")
    for exclusion in record["negative_results"]
        exclusion["status"] == "all_beta_rejected" || error("unresolved beta region")
        n = asint(exclusion["n0_cells"])
        first = seed_intervals(spectra, n, 1//big(1))
        inf(first.t) >= 1//2 || error("small-beta exclusion failed")
        tail = exact(exclusion["tail_beta"])
        en,rn = spectra[2n]; e2,r2 = spectra[4n]
        diff = minimum(Q(e)-r for (e,r) in zip(e2,r2)) -
            2minimum(Q(e)+r for (e,r) in zip(en,rn))
        diff > 0 || error("energy defect not positive")
        sup(exp(qinterval(4n)*log(interval(2))-qinterval(tail*diff))) <= 1//2 ||
            error("large-beta exclusion failed")
        previous = 1//big(1)
        for box in exclusion["cover"]
            a, b = exact(box["a"]), exact(box["b"])
            a == previous < b || error("beta cover has a gap or overlap")
            s, t = lower_defects(spectra[2n], spectra[4n], a, b)
            rejects(s,t) !== nothing || error("beta box is not excluded")
            previous = b
        end
        previous == tail || error("beta cover misses its tail")
    end
    return spectra
end

end
