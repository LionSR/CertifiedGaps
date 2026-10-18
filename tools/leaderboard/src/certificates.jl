# Reading the seed, defects and bound that a certificate file states, so that each leaderboard
# row can be compared with its source. The certificate formats differ by model but share a
# vocabulary of field names, listed in CERT_FIELDS. Every JSON object that names a seed length
# and an inverse temperature, together with a bound or exact defects, is a candidate; parameter
# fields of the object and its ancestors locate the parameter point. Objects with floating-point
# defects are diagnostic scans and are skipped.

using CertifiedGaps.Results: read_result

# Field names for each quantity, in order of preference (first present wins).
const CERT_FIELDS = (
    n0 = ("n0", "n0_cells", "n", "spatial_length", "even_length_threshold"),
    beta = ("beta0", "beta", "initial_beta"),
    bound = ("bound_exact_binary", "gap_lower_bound", "gap_lower", "gap_strict_lower", "bound"),
    s = ("epsilon_s", "eps_s", "eps_s_upper"),
    t = ("epsilon_t", "eps_t", "eps_t_upper"),
)
# Record parameter name => the certificate fields that may hold it. A Weyl-transported interval
# [D_min, D_max] is the certificate field D_interval.
const PARAM_FIELDS = Dict("g" => ("g",), "h" => ("h", "staggered_field"), "delta" => ("delta",),
    "Jperp" => ("J",), "D" => ("reference_D", "D"))
const SCOPE_FIELDS = union(Set(["D_interval"]), (Set(v) for v in values(PARAM_FIELDS))...)

"A JSON value as an exact rational: integers, exact strings, and floats converted exactly."
cert_exact(x::Union{AbstractFloat,Rational}) = Rational{BigInt}(x)
cert_exact(x) = parse_exact(x)

"Load a certificate (JSON, or gzipped JSON) as nested Dicts and Vectors."
load_certificate(path) = read_result(path, Dict{String,Any})

"""
    core_candidate(data) -> NamedTuple or nothing

The single entry of a certificate made by the shared pipeline (`GapCertificates.certify_gap`,
format "CertifiedGaps ring gap certificate v…"): the seed and exact defects from `seed`, and as
bound the outward-rounded interval bound `interval_bound_decimal` (`certified_gap_bound`), which
the model READMEs quote; the exact chain bound `gap_lower_bound` is at most slightly smaller.
"""
function core_candidate(data)
    format = get(data, "format", nothing)
    format isa AbstractString && startswith(format, "CertifiedGaps ring gap certificate") ||
        return nothing
    s = data["seed"]
    scope = Dict(k => v for (k, v) in data if k in SCOPE_FIELDS)
    return (; n0 = s["n0"], beta = s["beta0"], bound = data["interval_bound_decimal"],
        s = s["eps_s"], t = s["eps_t"], scope)
end

"""
    certificate_candidates(data) -> Vector{NamedTuple}

Every object in a certificate that states `n0`, `beta`, and a `bound` or exact defects, with
its defects `s`, `t` (or `nothing`) and the parameter fields in scope (its own and its
ancestors').
"""
function certificate_candidates(data)
    core = core_candidate(data)
    core === nothing || return [core]
    out = NamedTuple[]
    function walk(x, scope)
        if x isa AbstractVector
            foreach(y -> walk(y, scope), x)
        elseif x isa AbstractDict
            scope = merge(scope, Dict(k => v for (k, v) in x if k in SCOPE_FIELDS))
            pick(names) = (i = findfirst(k -> haskey(x, k), names); i === nothing ? nothing : x[names[i]])
            n0, beta, bound = pick(CERT_FIELDS.n0), pick(CERT_FIELDS.beta), pick(CERT_FIELDS.bound)
            # The spin-1 scalars list the rounded defects scale by scale; the seed is the last pair.
            rounded = get(x, "rounded_defects", nothing)
            defects = rounded isa AbstractVector && !isempty(rounded) ? rounded[end] : nothing
            bound isa Union{Real,AbstractString} || (bound = nothing)
            if n0 !== nothing && beta !== nothing && (bound !== nothing || defects isa AbstractVector)
                s = something(pick(CERT_FIELDS.s), defects isa AbstractVector ? defects[1] : nothing, Some(nothing))
                t = something(pick(CERT_FIELDS.t), defects isa AbstractVector ? defects[2] : nothing, Some(nothing))
                # Floating-point defects mark a diagnostic scan, not a certificate.
                (s isa AbstractFloat || t isa AbstractFloat) || push!(out, (; n0, beta, bound, s, t, scope))
            end
            foreach(v -> walk(v, scope), values(x))
        end
    end
    walk(data, Dict{String,Any}())
    return out
end

"Whether a candidate's parameter fields agree with a record's exact parameters."
function params_agree(scope, params)
    for (name, value) in params
        v = parse_exact(value)
        for field in get(PARAM_FIELDS, name, ())
            haskey(scope, field) && cert_exact(scope[field]) != v && return false
        end
        if name in ("D_min", "D_max") && haskey(scope, "D_interval")
            cert_exact(scope["D_interval"][name == "D_min" ? 1 : 2]) == v || return false
        end
    end
    return true
end

"""
    match_certificate(path, n0, beta0, params) -> candidate or nothing

The certificate entry for seed `(n0, beta0)` at the given parameter point, preferring the
largest stated bound when several entries qualify.
"""
function match_certificate(path, n0, beta0, params)
    best = nothing
    for c in certificate_candidates(load_certificate(path))
        (cert_exact(c.n0) == n0 && cert_exact(c.beta) == beta0 && params_agree(c.scope, params)) || continue
        better = best === nothing || (c.bound !== nothing && (best.bound === nothing || cert_exact(c.bound) > cert_exact(best.bound)))
        better && (best = c)
    end
    return best
end

"""
    check_certificate!(report, where, rec, headline, root)

Compare a record with the certificate it links: the certificate must state the record's seed
`(n0, beta0)` at its parameter point, and, unless the bound comes from a Lean theorem (a Lean
tier), a bound at least the record's (exact comparison).
Returns the matched certificate entry (for the Lean comparison) or `nothing`.
"""
function check_certificate!(report, where, rec, headline, root; path = rec["certificate"],
        seed = rec["seed"])
    full = joinpath(root, path)
    isfile(full) || return nothing
    n0, beta0 = parse_exact(seed["n0"]), parse_exact(seed["beta0"])
    c = try
        match_certificate(full, n0, beta0, something(get(rec, "params", nothing), Dict()))
    catch e
        err!(report, "$where: cannot read $path: $(sprint(showerror, e))")
        return nothing
    end
    if c === nothing
        err!(report, "$where: $path states no bound for seed ($(seed["n0"]), $(seed["beta0"])) at this parameter point")
    elseif c.bound === nothing
        startswith(rec["tier"], "lean-") ||
            err!(report, "$where: $path states no bound; only a Lean theorem can supply it")
    elseif headline !== nothing && !startswith(rec["tier"], "lean-") && cert_exact(c.bound) < headline[2]
# A Lean tier takes its bound from the Lean theorem, which may exceed the certificate's.
        err!(report, "$where: $path proves $(Float64(cert_exact(c.bound))), less than the row's bound")
    end
    return c
end

"""
    window_defects(path, n0, beta, params) -> (s, t) or nothing

The exact defects of a one-step entry of a rational chain certificate (the `finite_intervals` of
`CertifiedGaps.Recurrence.rational_gap_chain`, or the `bridge` of the spin-1 temperature
lifting): an object with `length_from = n0`, `beta = beta` and exact `p`, `q`, which bound
`1 - P_s(n0, β)` and `1 - P_τ(2n0, β)`. Errors if two such objects disagree.
"""
function window_defects(path, n0, beta, params)
    found = Set{Tuple{Rational{BigInt},Rational{BigInt}}}()
    function walk(x, scope)
        if x isa AbstractVector
            foreach(y -> walk(y, scope), x)
        elseif x isa AbstractDict
            scope = merge(scope, Dict(k => v for (k, v) in x if k in SCOPE_FIELDS))
            if all(k -> haskey(x, k), ("length_from", "beta", "p", "q")) &&
               x["p"] isa AbstractString && x["q"] isa AbstractString &&
               cert_exact(x["length_from"]) == n0 && cert_exact(x["beta"]) == beta &&
               params_agree(scope, params)
                push!(found, (parse_exact(x["p"]), parse_exact(x["q"])))
            end
            foreach(v -> walk(v, scope), values(x))
        end
    end
    walk(load_certificate(path), Dict{String,Any}())
    length(found) > 1 && throw(ArgumentError("$path: several defect pairs at ($n0, $beta)"))
    return isempty(found) ? nothing : first(found)
end
