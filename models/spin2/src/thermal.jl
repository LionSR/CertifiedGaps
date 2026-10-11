"""Approximate log Z_n(beta) from an MPO amplitude A = exp(-beta H/2).

Start from I/sqrt(5^n), evolve the ket legs, normalize each time step,
and accumulate log norms. Then log Z = n log 5 + 2 log(norm growth).
A symmetric forward/reverse bond sweep has second-order Trotter error.
No ground-state projection or restriction of the thermal trace to Sz=0.
Implemented by `CertifiedGaps.Purification.purified_log_partition`.
"""
function thermal_curve(
    n,
    betas;
    dbeta = 0.05,
    maxdim = 64,
    cutoff = 1e-10,
    periodic = true,
    progress = true,
)
    pairs = bonds(n; periodic)
    dbeta > 0 || throw(ArgumentError("dbeta must be positive"))
    maxdim > 0 || throw(ArgumentError("maxdim must be positive"))
    sites = siteinds("S=2", n; conserve_qns = true)
    halfgate = (sites, j, k, tau) -> exp(-tau*bond_operator(sites, j, k))
    return Purification.purified_log_partition(
        sites,
        pairs,
        halfgate,
        betas;
        dbeta,
        maxdim,
        cutoff,
        progress,
    )
end

"""Evaluate both seed ratios without clamping apparent violations to [0,1]."""
function seed_estimates(short, long, n, beta)
    logs = long[beta]["logZ"]-2short[beta]["logZ"]
    logt = long[2beta]["logZ"]-2long[beta]["logZ"]
    ps, pt = exp(logs), exp(logt)
    return Dict(
        "n"=>n,
        "beta"=>beta,
        "spatial_purity"=>ps,
        "thermal_purity_2n"=>pt,
        "log_spatial_purity"=>logs,
        "log_thermal_purity_2n"=>logt,
        "numerical_seed_candidate"=>0.99 <= ps <= 1 && 0.99 <= pt <= 1,
        "physical_range_ok"=>0 < ps <= 1 && 0 < pt <= 1,
        "certified"=>false,
    )
end
