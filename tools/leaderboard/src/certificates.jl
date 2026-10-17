# Reading the seed, defects and bound that a certificate file states, so that each leaderboard
# row can be compared with its source. The certificate formats differ by model but share a
# vocabulary of field names, listed in CERT_FIELDS. Every JSON object that names a seed length
# and an inverse temperature, together with a bound or exact defects, is a candidate; parameter
# fields of the object and its ancestors locate the parameter point. Objects with floating-point
# defects are diagnostic scans and are skipped.

using CodecZlib: GzipDecompressor, transcode

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
const PARAM_FIELDS = Dict("g" => ("g",), "h" => ("h",), "delta" => ("delta",), "Jperp" => ("J",),
    "D" => ("reference_D", "D"))
const SCOPE_FIELDS = union(Set(["D_interval"]), (Set(v) for v in values(PARAM_FIELDS))...)

"A JSON value as an exact rational: integers, exact strings, and floats converted exactly."
cert_exact(x::AbstractFloat) = Rational{BigInt}(x)
cert_exact(x) = parse_exact(x)

"Load a certificate (JSON, or gzipped JSON) as nested Dicts and Vectors."
function load_certificate(path)
    bytes = read(path)
    endswith(path, ".gz") && (bytes = transcode(GzipDecompressor, bytes))
    return JSON3.read(bytes, Dict{String,Any})
end

"""
    certificate_candidates(data) -> Vector{NamedTuple}

Every object in a certificate that states `n0`, `beta`, and a `bound` or exact defects, with
its defects `s`, `t` (or `nothing`) and the parameter fields in scope (its own and its
ancestors').
"""
function certificate_candidates(data)
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
`(n0, beta0)` at its parameter point, and a bound at least the record's (exact comparison).
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
    elseif headline !== nothing && cert_exact(c.bound) < headline[2]
        err!(report, "$where: $path proves $(Float64(cert_exact(c.bound))), less than the row's bound")
    end
    return c
end
