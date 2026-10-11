"""
    ThermalMPOGap.Purification

Numerical thermal purification: log Z(β) of a spin system from an MPO amplitude
A = exp(−βH/2), evolved from the normalised identity with a second-order product formula.
This is a numerical control, never a certificate.

The method is implemented by the ITensors extension and becomes available after
`using ITensors, ITensorMPS`.
"""
module Purification

export purified_log_partition

"""
    purified_log_partition(sites, pairs, halfgate, betas; dbeta, maxdim, cutoff,
        apply_kwargs = (;), extra = nothing, progress = false, label = "", callback = nothing)

Approximate log Z(β) for every β in `betas` (each a multiple of `dbeta`).

* `sites`: ITensors site indices; `pairs`: the bonds `(i, j)`;
* `halfgate(sites, i, j, τ)`: the ITensor `exp(−τ h_ij)`; each step applies the half-gates for
  τ = dbeta/4 forward, then in reverse (a symmetric second-order step of size dbeta/2 on A);
* `apply_kwargs`: extra keywords for `ITensors.apply` (e.g. `move_sites_back_between_gates`);
* `extra(logZ)`: additional fields for each row with β > 0.

Starts from I/√(dⁿ), renormalises after every step and returns
`Dict(β => Dict("beta", "logZ", "maxlinkdim", "seconds", extra...))` with
log Z = n log d + 2 Σ log(norm growth). `callback(result)` runs after each target β.
"""
function purified_log_partition end

end
