"""
    CertifiedGaps.Exclusions

Negative results for the doubling recurrence: proofs that a seed length n fails at every β > 0,
from the verified spectra of the rings n and 2n. Lower bounds on the seed defects are pushed
through the order-preserving doubling map until an argument z reaches 1/2, after which both
components stay ≥ 1/2 and the recurrence cannot converge to zero (`rejection`); the β axis is
covered by bisected windows [lo, hi] with such a witness and a tail [B, ∞) bounded through the
ground-energy mismatch (`exclude_all_beta`). A rejection is never a gaplessness claim: it rules
out this seed length for the stated recurrence only.
"""
module Exclusions

using IntervalArithmetic
using ..ExactArithmetic: RationalBig, round_down, parse_rational
using ..Recurrence: to_interval
using ..VerifiedSpectra: BlockSpectrum, shifted_log_partition, ground_energy_enclosure

export rejection, exclude_all_beta, check_exclusion

const Q = RationalBig
const GRID = big(10)^15                  # lower defects are rounded down to this grid
const SCOPE = "No positive rate from the doubling recurrence at this seed length for any " *
              "beta > 0; not a gaplessness claim"

_token(x) = (x = Q(x); denominator(x) == 1 ? string(numerator(x)) :
                       string(numerator(x), "//", denominator(x)))

"""
    rejection(s, t; steps = 64) -> Dict or nothing

A witness that lower bounds `(s, t)` ≤ (1 − Pₛ, 1 − P_τ) on the seed defects rule out every
gap bound from the doubling recurrence: iterate the map on lower bounds (rounded down, which the
order-preserving map tolerates) until an argument z reaches 1/2, after which both components
stay ≥ 1/2. `nothing` means no witness was found, never a certificate.
"""
function rejection(s::Q, t::Q; steps = 64)
    0 <= s <= 1 && 0 <= t <= 1 || throw(ArgumentError("invalid lower defects"))
    grid = big(10)^24
    for j in 0:steps
        z = 1 - (1 - t) * (1 - s)^2
        # z ≥ 1/2 implies both next components ≥ 1/2, hence no attraction to zero.
        z >= 1 // 2 && return Dict("step" => _token(j), "z_lower" => _token(z))
        sp = round_down(z^2 / (2 * (1 - z)^2), grid)
        zt = 1 - (1 - sp) * (1 - t)^2
        zt >= 1 // 2 && return Dict("step" => _token(j), "thermal_z_lower" => _token(zt))
        s, t = sp, round_down(zt^2 / (2 * (1 - zt)^2), grid)
    end
    return nothing
end
rejection(s, t; kwargs...) = rejection(Q(s), Q(t); kwargs...)

_lower(x::Interval) = max(zero(Q), round_down(Q(inf(x)), GRID))

# Exact lower ends a, b of the ground energies E₀(n), E₀(2n) (as tuples (lo, hi)), and a cached
# F_L(β) = log Tr exp[−β(H_L − a_L)] with a_L the lower end of E₀(L), so that H_L − a_L ⪰ 0 and
# F_L decreases in β.
function _data(small, large)
    a, b = ground_energy_enclosure(small), ground_energy_enclosure(large)
    cache = Dict{Tuple{Bool,Q},Interval{Float64}}()
    F(isbig, β) = get!(cache, (isbig, Q(β))) do
        isbig ? shifted_log_partition(large, Q(β), b[1]) : shifted_log_partition(small, Q(β), a[1])
    end
    return (; a, b, F)
end

# Lower defects valid for every β in [lo, hi] (Float64-exact endpoints). Spatial: F_L
# decreases, so log Pₛ ≤ F_2n(lo) − 2F_n(hi) − β(a_2n − 2a_n) over the whole interval; no
# monotonicity of Pₛ is assumed. Thermal: d log P_τ/dβ = 2(⟨H⟩_β − ⟨H⟩_2β) ≥ 0, so
# P_τ(2n, β) ≤ P_τ(2n, hi).
function _window(data, lo::Q, hi::Q)
    Q(Float64(lo)) == lo && Q(Float64(hi)) == hi || throw(ArgumentError("non-dyadic window"))
    β = interval(Float64(lo), Float64(hi))
    logps = data.F(true, lo) - 2data.F(false, hi) - β * to_interval(data.b[1] - 2data.a[1])
    logpt = data.F(true, 2hi) - 2data.F(true, hi)
    return _lower(-expm1(logps)), _lower(-expm1(logpt))
end

# Rejection of every β ≥ B: with δ = lower E₀(2n) − 2·upper E₀(n) > 0,
# log Pₛ(β) ≤ F_2n(β) − βδ ≤ F_2n(B) − Bδ.
function _tail(data, B::Q)
    δ = data.b[1] - 2data.a[2]
    δ > 0 || return nothing
    s = _lower(-expm1(data.F(true, B) - to_interval(B * δ)))
    rejection(s, zero(Q)) === nothing && return nothing
    return Dict("beta_from" => _token(B), "ground_mismatch_lower" => _token(δ),
        "eps_s_lower" => _token(s), "eps_t_lower" => "0")
end

"""
    exclude_all_beta(small, large, n; max_depth = 20, max_tail = 1024) -> Dict

Proof that the seed length `n` fails the doubling recurrence at every β > 0, from the verified
spectra `small`, `large` of the rings n and 2n (`verified_ring_spectrum`): a tail start B (the
smallest power of 2 up to `max_tail` whose tail is rejected), then [0, B] bisected into windows
that each carry a `rejection` witness of their lower defects. Throws if a window stays
unresolved after `max_depth` bisections.
"""
function exclude_all_beta(small::AbstractVector{BlockSpectrum},
    large::AbstractVector{BlockSpectrum}, n::Integer; max_depth = 20, max_tail = 1024)
    data = _data(small, large)
    B = one(Q)
    tail = _tail(data, B)
    while tail === nothing && B < max_tail
        B *= 2
        tail = _tail(data, B)
    end
    tail === nothing && error("n = $n: the large-β tail is not rejected up to β = $max_tail")
    windows = Any[]
    function split(lo, hi, depth)
        s, t = _window(data, lo, hi)
        witness = rejection(s, t)
        if witness !== nothing
            push!(windows, Dict("beta_from" => _token(lo), "beta_to" => _token(hi),
                "eps_s_lower" => _token(s), "eps_t_lower" => _token(t), "rejection" => witness))
        else
            depth < max_depth || error("n = $n: β window [$lo, $hi] unresolved")
            split(lo, (lo + hi) / 2, depth + 1)
            split((lo + hi) / 2, hi, depth + 1)
        end
    end
    split(zero(Q), B, 0)
    return Dict{String,Any}("n0" => _token(n), "windows" => windows, "tail" => tail,
        "scope" => SCOPE)
end

"""
    check_exclusion(record; small = nothing, large = nothing) -> true

Check an `exclude_all_beta` record: the windows tile [0, B] without gaps, the tail starts at B,
and every stored lower defect carries its rejection witness. With the verified spectra `small`,
`large` of the rings n and 2n, every lower defect and the ground-energy mismatch are also
recomputed and the stored values must not exceed them. Throws on any failure.
"""
function check_exclusion(record; small = nothing, large = nothing)
    data = small === nothing ? nothing : _data(small, large)
    q(x) = parse_rational(string(x))
    previous = zero(Q)
    for w in record["windows"]
        lo, hi = q(w["beta_from"]), q(w["beta_to"])
        lo == previous && hi > lo || error("β coverage gap or overlap at $lo")
        s, t = q(w["eps_s_lower"]), q(w["eps_t_lower"])
        rejection(s, t) == w["rejection"] || error("rejection witness fails on [$lo, $hi]")
        if data !== nothing
            fs, ft = _window(data, lo, hi)
            s <= fs && t <= ft || error("stored lower defects exceed the fresh ones on [$lo, $hi]")
        end
        previous = hi
    end
    tail = record["tail"]
    q(tail["beta_from"]) == previous || error("the tail does not start where the windows end")
    q(tail["eps_t_lower"]) == 0 || error("invalid tail thermal defect")
    rejection(q(tail["eps_s_lower"]), zero(Q)) !== nothing || error("tail witness fails")
    if data !== nothing
        fresh = _tail(data, previous)
        fresh !== nothing &&
            q(tail["ground_mismatch_lower"]) <= q(fresh["ground_mismatch_lower"]) &&
            q(tail["eps_s_lower"]) <= q(fresh["eps_s_lower"]) || error("the tail bound fails")
    end
    return true
end

end
