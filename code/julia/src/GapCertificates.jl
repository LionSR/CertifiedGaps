"""
    CertifiedGaps.GapCertificates

End-to-end gap certificates for a `RingModel` whose seed rings are small enough for exact
diagonalisation: the exact RP identity, rigorous seed defects from verified spectra of the
rings n₀ and 2n₀, and the exact rational doubling chain. The certificate is a JSON-ready
dictionary in which every exact number is a string; `check_certificate` re-derives it.

The seed is stored as (n₀, β₀, εₛ, εₜ), all exact, with εₛ ≥ 1 − Pₛ(n₀, β₀) and
εₜ ≥ 1 − P_τ(2n₀, β₀), and the RP data as −h = Σ w_α C_α ⊗ conj(C_α): the hypotheses of the
Lean ring theorem.
"""
module GapCertificates

using IntervalArithmetic
using ..ExactArithmetic: RationalBig, round_up, parse_rational
using ..RingModels
using ..RingModels: _int
using ..VerifiedSpectra
using ..Recurrence: rational_gap_chain, certified_gap_bound
using ..Seeds: seed_defects
using ..Results: read_result, write_result, model_metadata, verify_model_metadata

export seed_bounds, certify_gap, check_certificate, certify_point, check_point,
    certificate_row, lean_points

const FORMAT = "CertifiedGaps ring gap certificate v1"
const GRID = big(10)^15                          # seed defects are rounded up to this grid

_upper(x::Interval) = max(round_up(RationalBig(sup(x)), GRID), 1 // GRID)
_endpoints(x::Interval) = [string(RationalBig(inf(x))), string(RationalBig(sup(x)))]

# The model whose rings are diagonalized and its cell size: the model itself, or the exact cell
# model of a `blocking = (; cell, charge, kwargs...)` (keywords of `cell_model`).
_spectral(model, ::Nothing) = model, 1
_spectral(model, b::NamedTuple) =
    cell_model(model, b.cell; Base.structdiff(b, NamedTuple{(:cell,)})...), b.cell

_blocking_record(::Nothing, spectral, c) = nothing
_blocking_record(::NamedTuple, spectral, c) =
    Dict{String,Any}("cell" => c, "model" => encode_model(spectral))

"""
    seed_bounds(model, n0, β0; blocking = nothing, positive_transfer = "", spectra = nothing,
        progress = false) -> NamedTuple

Rigorous seed defects at `n0` and exact `β0`: enclosures of log Z_{n₀}(β₀), log Z_{2n₀}(β₀)
and log Z_{2n₀}(2β₀) from `verified_ring_spectrum`, and the defect enclosures

    1 − Pₛ(n₀, β₀) = 1 − exp(log Z_{2n₀}(β₀) − 2 log Z_{n₀}(β₀)),
    1 − P_τ(2n₀, β₀) = 1 − exp(log Z_{2n₀}(2β₀) − 2 log Z_{2n₀}(β₀)),

with `eps_s`, `eps_t` their upper ends rounded up to a rational grid. The verified spectra of
the two rings are returned as `small` and `large`.

* `blocking = (; cell = c, charge, modulus, flip, …)` diagonalizes the rings as rings of
  n₀/c and 2n₀/c cells of `cell_model(model, c; charge, …)`: the same Hamiltonians, split by a
  charge (or symmetries) that only the cells carry.
* An odd `n0` needs `positive_transfer`: a reference to the model's proof that its transfer
  operator is positive, which the criterion then uses in place of an even seed length.
* `spectra = Dict(n0 => small, 2n0 => large)` (keyed by ring length in sites) reuses spectra
  from an earlier `seed_bounds` or `verified_ring_spectrum` call of the same (cell) model, for
  instance over a grid of β₀.
"""
function seed_bounds(model::RingModel, n0::Integer, β0; blocking = nothing,
    positive_transfer::AbstractString = "", spectra = nothing, progress = false)
    n0 >= 2 || throw(ArgumentError("n0 must be at least 2"))
    iseven(n0) || !isempty(positive_transfer) ||
        throw(ArgumentError("an odd n0 needs a positive_transfer reference"))
    β = RationalBig(β0)
    β > 0 || throw(ArgumentError("β0 must be positive"))
    spectral, c = _spectral(model, blocking)
    n0 % c == 0 || throw(ArgumentError("n0 = $n0 is not a whole number of $c-site cells"))
    ring(L) = spectra === nothing ? verified_ring_spectrum(spectral, L ÷ c; progress) : spectra[L]
    small, large = ring(n0), ring(2n0)
    zn = log_partition_enclosure(small, β)
    z2n = log_partition_enclosure(large, β)
    z2n2 = log_partition_enclosure(large, 2β)
    (; defect_s, defect_t) = seed_defects(zn, z2n, z2n2)
    return (; n0, β0 = β, blocking = _blocking_record(blocking, spectral, c), positive_transfer,
        eps_s = _upper(defect_s), eps_t = _upper(defect_t), defect_s, defect_t,
        log_partition = Dict("logZ_n0_beta0" => _endpoints(zn),
            "logZ_2n0_beta0" => _endpoints(z2n), "logZ_2n0_2beta0" => _endpoints(z2n2)),
        spectra = Dict(string(n0) => spectrum_summary(small),
            string(2n0) => spectrum_summary(large)),
        small, large)
end

# The exact chain with the number of doublings that gives the best bound. Each defect is rounded
# up to the grid 10⁻⁵¹², so once the true defects fall below it the tail bound log(1/(C u))/β_J
# is limited by the grid rather than the gap; too few doublings and the tail is not yet reached.
function _chain(εs, εt, n0, β; max_steps = 16)
    best = nothing
    for steps in 1:max_steps
        chain = try
            rational_gap_chain(εs, εt, n0, β; steps)
        catch err
            err isa AssertionError || rethrow()
            continue
        end
        if best === nothing || parse_rational(chain["gap_strict_lower"]) >
                               parse_rational(best["gap_strict_lower"])
            best = chain
        end
    end
    best === nothing &&
        error("the rational chain does not close for εₛ = $(Float64(εs)), εₜ = $(Float64(εt))")
    return best
end

"""
    certify_gap(model; n0, β0, rp = rp_certificate(model), blocking = nothing,
        positive_transfer = "", seed = nothing, params = nothing, progress = false) -> Dict

Gap certificate for every even ring length L ≥ `n0`: verifies the RP certificate `rp`, bounds
the seed defects with `seed_bounds` (with `blocking` and `positive_transfer` as there, both
recorded in the certificate), and runs `rational_gap_chain`. A `seed` already computed by
`seed_bounds` with the same arguments is used instead of recomputing the spectra. Throws if the
seed is not admissible (the interval bound `certified_gap_bound` is the screen). The result's
`"gap_lower_bound"` is a strict lower bound on Δ_L. `params`, the exact parameters of the
model's point (a `Dict` or `NamedTuple` of integers and rationals), are recorded as `"params"`.
"""
function certify_gap(model::RingModel; n0::Integer, β0, rp = rp_certificate(model),
    blocking = nothing, positive_transfer::AbstractString = "", seed = nothing,
    params = nothing, progress = false)
    verify_rp(model, rp) || error("the RP certificate does not verify")
    β = RationalBig(β0)
    if seed === nothing
        seed = seed_bounds(model, n0, β; blocking, positive_transfer, progress)
    else
        spectral, c = _spectral(model, blocking)
        (seed.n0, seed.β0, seed.positive_transfer) == (n0, β, positive_transfer) &&
            seed.blocking == _blocking_record(blocking, spectral, c) ||
            throw(ArgumentError("the given seed was computed for other arguments"))
    end
    screen = certified_gap_bound(seed.eps_s, seed.eps_t, β)
    screen.bound === nothing && error("seed not admissible: εₛ ≤ $(Float64(seed.eps_s)), " *
                                      "εₜ ≤ $(Float64(seed.eps_t)) at (n₀, β₀) = ($n0, $β)")
    chain = _chain(seed.eps_s, seed.eps_t, n0, β)
    record = Dict{String,Any}("n0" => n0, "beta0" => string(β), "eps_s" => string(seed.eps_s),
        "eps_t" => string(seed.eps_t))
    isempty(positive_transfer) || (record["positive_transfer"] = positive_transfer)
    cert = Dict{String,Any}(
        "format" => FORMAT,
        "claim" => "every even ring length L >= n0 has a unique ground state and " *
                   "spectral gap > gap_lower_bound",
        "model" => encode_model(model),
        "reflection_positivity" => encode_rp(rp),
        "seed" => record,
        "log_partition" => seed.log_partition,
        "spectra" => seed.spectra,
        "chain" => chain,
        "gap_lower_bound" => chain["gap_strict_lower"],
        "gap_lower_bound_decimal" => chain["gap_strict_lower_decimal"],
        "interval_bound_decimal" => screen.bound,
    )
    seed.blocking === nothing || (cert["blocking"] = seed.blocking)
    params === nothing || (cert["params"] = _params(params))
    return cert
end

# Exact parameters as strings "num//den" (or "int"); floats are refused.
_params(params) = Dict{String,Any}(string(k) => _exact_string(v) for (k, v) in pairs(params))
_exact_string(x::Union{Integer,Rational}) = denominator(x) == 1 ? string(numerator(x)) :
    string(numerator(x), "//", denominator(x))
_exact_string(x::AbstractString) = _exact_string(parse_rational(x))
_exact_string(x) = throw(ArgumentError("refusing inexact parameter $(repr(x))"))

# The `blocking` argument that reproduces a certificate's stored cell model from its site model.
function _stored_blocking(cert, model)
    haskey(cert, "blocking") || return nothing
    c = _int(cert["blocking"]["cell"])
    stored = decode_model(cert["blocking"]["model"])
    blocking = (; cell = c, charge = stored.charge, modulus = stored.modulus, flip = stored.flip,
        symmetries = stored.symmetries, reflection = stored.reflection)
    rebuilt, _ = _spectral(model, blocking)
    rebuilt.bond == stored.bond && rebuilt.onsite == stored.onsite ||
        error("the stored cell model is not the blocked site model")
    return blocking
end

"""
    check_certificate(cert; spectra = true, progress = false) -> Dict

Re-derive a certificate from its stored model, RP data and seed: the exact RP identity, the
seed defects from freshly verified spectra (skipped with `spectra = false`), which must not
exceed the stored εₛ, εₜ, and the exact chain, whose bound must equal the stored one. Throws
on any mismatch; returns the recomputed defect enclosures. A checker that has itself just
computed the verified spectra of the certificate's (cell) model, from the decoded model, may
pass them as `spectra = Dict(n0 => small, 2n0 => large)` (keyed by length in sites) instead of
having them recomputed; spectra read from storage are not a check.
"""
function check_certificate(cert; spectra = true, progress = false)
    cert["format"] == FORMAT || error("unknown certificate format $(cert["format"])")
    model = decode_model(cert["model"])
    rp = decode_rp(cert["reflection_positivity"])
    verify_rp(model, rp) || error("RP identity fails")
    s = cert["seed"]
    n0 = _int(s["n0"])
    β, εs, εt = parse_rational(s["beta0"]), parse_rational(s["eps_s"]), parse_rational(s["eps_t"])
    positive_transfer = String(get(s, "positive_transfer", ""))
    iseven(n0) || !isempty(positive_transfer) || error("odd n0 without positive_transfer")
    blocking = _stored_blocking(cert, model)
    recomputed = Dict{String,Any}()
    if spectra !== false
        fresh = seed_bounds(model, n0, β; blocking, positive_transfer, progress,
            spectra = spectra === true ? nothing : spectra)
        sup(fresh.defect_s) <= εs || error("spatial defect exceeds the stored εₛ")
        sup(fresh.defect_t) <= εt || error("thermal defect exceeds the stored εₜ")
        recomputed["defect_s_upper"] = sup(fresh.defect_s)
        recomputed["defect_t_upper"] = sup(fresh.defect_t)
    end
    chain = rational_gap_chain(εs, εt, n0, β; steps = length(cert["chain"]["finite_intervals"]))
    chain["gap_strict_lower"] == cert["gap_lower_bound"] == cert["chain"]["gap_strict_lower"] ||
        error("the rational chain does not reproduce the stored bound")
    recomputed["gap_lower_bound"] = chain["gap_strict_lower"]
    screen = certified_gap_bound(εs, εt, β).bound
    screen !== nothing && screen == _float(cert["interval_bound_decimal"]) ||
        error("the interval bound does not reproduce the stored interval_bound_decimal")
    recomputed["interval_bound_decimal"] = screen
    return recomputed
end

_float(x::AbstractString) = parse(Float64, x)
_float(x::Real) = Float64(x)

"""
    certify_point(path, model; n0, β0, params, model_dir, sources = (;), packages = (),
        record = Dict(), kwargs...) -> Dict

The certificate of one parameter point of a model, written to `path`: `certify_gap(model; n0,
β0, params, kwargs...)` (`rp`, `blocking`, `positive_transfer`, `seed`, `progress`), the extra
top-level fields `record`, and under `"metadata"` the run record `Results.model_metadata` of
`model_dir` with the model-local `sources` (`dirs`, `files`), the configuration (the params,
n₀, β₀ and the run time) and the versions of `CertifiedGaps` and `packages`.
"""
function certify_point(path, model::RingModel; n0::Integer, β0, params, model_dir,
    sources = (;), packages = (), record = Dict{String,Any}(), kwargs...)
    seconds = @elapsed cert = certify_gap(model; n0, β0, params, kwargs...)
    merge!(cert, record)
    config = merge(_params(params), Dict{String,Any}("n0" => n0,
        "beta0" => _exact_string(RationalBig(β0)), "seconds" => seconds))
    cert["metadata"] = model_metadata(model_dir, config; sources...,
        packages = (parentmodule(@__MODULE__), packages...), certified = true)
    write_result(path, cert)
    return cert
end

"""
    check_point(path, model; n0, β0, params, model_dir, sources = (;), rp = rp_certificate(model),
        blocking = nothing, positive_transfer = "", lean = nothing, spectra = false) -> Dict

Check the stored certificate `path` of one parameter point against its definition, then run
`check_certificate(cert; spectra)`. Requires the recorded model-local source hashes
(`Results.verify_model_metadata` of `model_dir` with `sources`), and exactly the point's model
record (`encode_model(model)`), RP record (`encode_rp(rp)`), seed (n₀, β₀, `positive_transfer`),
cell blocking and `params`. `lean = (s, t)` are the seed defects hardcoded in the model's Lean
theorem (see `lean_points`), which the stored εₛ, εₜ must not exceed. Throws on any mismatch
and returns the values recomputed by `check_certificate`.
"""
function check_point(path, model::RingModel; n0::Integer, β0, params, model_dir,
    sources = (;), rp = rp_certificate(model), blocking = nothing,
    positive_transfer::AbstractString = "", lean = nothing, spectra = false)
    name = basename(path)
    fail(what) = error("$name: $what")
    cert = read_result(path)
    verify_model_metadata(cert["metadata"], model_dir; sources...)
    encode_model(decode_model(cert["model"])) == encode_model(model) ||
        fail("the stored model is not the point's model")
    encode_rp(decode_rp(cert["reflection_positivity"])) == encode_rp(rp) ||
        fail("the stored RP data are not the point's")
    s = cert["seed"]
    (_int(s["n0"]), parse_rational(s["beta0"])) == (n0, RationalBig(β0)) ||
        fail("the stored seed (n0, beta0) is not the point's")
    String(get(s, "positive_transfer", "")) == positive_transfer ||
        fail("unexpected positive_transfer reference")
    spectral, c = _spectral(model, blocking)
    stored = haskey(cert, "blocking") ?
        (_int(cert["blocking"]["cell"]), encode_model(decode_model(cert["blocking"]["model"]))) :
        nothing
    stored == (blocking === nothing ? nothing : (c, encode_model(spectral))) ||
        fail("the stored blocking is not the point's")
    Dict(String(k) => _exact_string(v) for (k, v) in get(cert, "params", Dict())) ==
        _params(params) || fail("the stored params are not the point's")
    if lean !== nothing
        parse_rational(s["eps_s"]) <= lean.s && parse_rational(s["eps_t"]) <= lean.t ||
            fail("the seed defects exceed those of the Lean theorem")
    end
    return check_certificate(cert; spectra)
end

"""
    certificate_row(cert) -> Dict

Summary row of a certificate: its params, the seed (n₀, β₀, εₛ, εₜ) and the bounds.
"""
function certificate_row(cert)
    row = Dict{String,Any}(String(k) => v for (k, v) in get(cert, "params", Dict()))
    for key in ("n0", "beta0", "eps_s", "eps_t")
        row[key] = cert["seed"][key]
    end
    for key in ("gap_lower_bound", "gap_lower_bound_decimal", "interval_bound_decimal")
        row[key] = cert[key]
    end
    return row
end

"""
    lean_points(path) -> Vector{NamedTuple}

The points of a Lean model file (`lean/CertifiedGaps/Models/*.lean`): each entry
`⟨p₁, …, β₀, { s := [s₀, …] t := [t₀, …] … }⟩` of its `points` list as `(; params, s, t)` with
the exact numbers `params = [p₁, …, β₀]` and the seed defects `s₀`, `t₀` that the Lean theorem
takes as hypotheses. Numbers are Lean ℚ literals: `7`, `13 / 2` or `9064e-9`.
"""
function lean_points(path)
    text = read(path, String)
    start = findfirst("def points", text)
    start === nothing && throw(ArgumentError("$path has no `def points` list"))
    number(x) = (x = replace(strip(x), " " => "");
        m = match(r"^(\d+)e-(\d+)$", x);
        m !== nothing ? parse(BigInt, m[1]) // big(10)^parse(Int, m[2]) : parse_rational(x))
    entry = r"⟨([^{⟩]*),\s*\{\s*s\s*:=\s*\[([^,\]]+)[^\]]*\]\s*t\s*:=\s*\[([^,\]]+)"
    points = [(; params = number.(split(m[1], ",")), s = number(m[2]), t = number(m[3]))
              for m in eachmatch(entry, text[last(start):end])]
    isempty(points) && throw(ArgumentError("no points found in $path"))
    return points
end

end
