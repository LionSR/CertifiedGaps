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
using ..VerifiedSpectra
using ..Recurrence: rational_gap_chain, certified_gap_bound

export seed_bounds, certify_gap, check_certificate

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
    defect_s = 1 - exp(z2n - 2zn)
    defect_t = 1 - exp(z2n2 - 2z2n)
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
        positive_transfer = "", seed = nothing, progress = false) -> Dict

Gap certificate for every even ring length L ≥ `n0`: verifies the RP certificate `rp`, bounds
the seed defects with `seed_bounds` (with `blocking` and `positive_transfer` as there, both
recorded in the certificate), and runs `rational_gap_chain`. A `seed` already computed by
`seed_bounds` with the same arguments is used instead of recomputing the spectra. Throws if the
seed is not admissible (the interval bound `certified_gap_bound` is the screen). The result's
`"gap_lower_bound"` is a strict lower bound on Δ_L.
"""
function certify_gap(model::RingModel; n0::Integer, β0, rp = rp_certificate(model),
    blocking = nothing, positive_transfer::AbstractString = "", seed = nothing,
    progress = false)
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
    return cert
end

# The `blocking` argument that reproduces a certificate's stored cell model from its site model.
function _stored_blocking(cert, model)
    haskey(cert, "blocking") || return nothing
    c = _int(cert["blocking"]["cell"])
    stored = decode_model(cert["blocking"]["model"])
    blocking = (; cell = c, charge = stored.charge, modulus = stored.modulus, flip = stored.flip)
    rebuilt, _ = _spectral(model, blocking)
    rebuilt.bond == stored.bond && rebuilt.onsite == stored.onsite ||
        error("the stored cell model is not the blocked site model")
    return blocking
end
_int(x::AbstractString) = parse(Int, x)
_int(x::Integer) = Int(x)

"""
    check_certificate(cert; spectra = true, progress = false) -> Dict

Re-derive a certificate from its stored model, RP data and seed: the exact RP identity, the
seed defects from freshly verified spectra (skipped with `spectra = false`), which must not
exceed the stored εₛ, εₜ, and the exact chain, whose bound must equal the stored one. Throws
on any mismatch; returns the recomputed defect enclosures.
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
    if spectra
        fresh = seed_bounds(model, n0, β; blocking, positive_transfer, progress)
        sup(fresh.defect_s) <= εs || error("spatial defect exceeds the stored εₛ")
        sup(fresh.defect_t) <= εt || error("thermal defect exceeds the stored εₜ")
        recomputed["defect_s_upper"] = sup(fresh.defect_s)
        recomputed["defect_t_upper"] = sup(fresh.defect_t)
    end
    chain = rational_gap_chain(εs, εt, n0, β; steps = length(cert["chain"]["finite_intervals"]))
    chain["gap_strict_lower"] == cert["gap_lower_bound"] == cert["chain"]["gap_strict_lower"] ||
        error("the rational chain does not reproduce the stored bound")
    recomputed["gap_lower_bound"] = chain["gap_strict_lower"]
    return recomputed
end

end
