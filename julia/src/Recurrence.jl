"""
    ThermalMPOGap.Recurrence

The two-purity doubling criterion. Seed defects εₛ = 1 − Pₛ(n₀, β₀) and
εₜ = 1 − P_τ(2n₀, β₀) are pushed through the exact doubling map; the per-scale bounds
log((1 − z_j)/z_j)/β_j give a uniform gap for every even length ≥ n₀.

Three evaluations of the same map, by number type:

* `map.jl`: any real type (`Float64`, `BigFloat`, intervals); the threshold ε⋆ and the
  majorant constants of survey §2.6;
* `certified.jl`: outward-rounded interval arithmetic on log-defects, with the r_J ≤ 1/10 tail
  (survey Lemma 2.21);
* `rational.jl`: exact rationals rounded up to a grid, with the bootstrap tail of the source
  paper (as in the spin-1 certificates).
"""
module Recurrence

using IntervalArithmetic
using ..ExactArithmetic: RationalBig, frac, round_up, round_down, log_enclosure

include("recurrence/map.jl")
include("recurrence/certified.jl")
include("recurrence/rational.jl")

end
