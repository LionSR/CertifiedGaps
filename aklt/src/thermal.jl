"""All-sector purification A=e^(-beta H/2) using exact projector gates.

Only the physical ket is evolved. The ancilla/bra indices retain the full trace.
Symmetric forward/reverse bond sweeps give a second-order product formula.
Norms reconstruct log Z, without resetting the ground energy or clipping R.
The evolution is the core `purified_log_partition`; each row also records
`R = Z-1` and whether `log Z >= 0`.
"""
function thermal_curve(
    n,
    bonds,
    betas;
    dbeta = 0.125,
    maxdim = 32,
    cutoff = 1e-12,
    move_sites_back_between_gates = true,
    progress = true,
    callback = nothing,
)
    sites = siteinds("AKLT3/2", n; conserve_qns = true)
    return purified_log_partition(
        sites,
        bonds,
        bond_gate,
        betas;
        dbeta,
        maxdim,
        cutoff,
        apply_kwargs = (; move_sites_back_between_gates),
        extra = logz->Dict("R"=>expm1(logz), "zero_energy_range_ok"=>logz>=0),
        progress,
        callback,
    )
end

function log_excess(logz)
    logz>0 || return nothing
    return logz>log(2) ? logz+log1p(-exp(-logz)) : log(expm1(logz))
end

"""Numerical diagnostics only; positive gluing defects are preserved."""
function diagnostics(small, wide, large, beta)
    ls = small[beta]["logZ"]
    ls2 = small[2beta]["logZ"]
    lw2 = wide[2beta]["logZ"]
    ll2 = large[2beta]["logZ"]
    lr, lrnext = log_excess(ls), log_excess(ll2)
    logc = isnothing(lr) || isnothing(lrnext) ? nothing : lrnext-2lr
    g = ll2-4ls2
    eta = max(0.0, g)
    eta_over_r2 = isnothing(lr) ? nothing : eta*exp(-2lr)
    return Dict(
        "beta"=>beta,
        "R"=>expm1(ls),
        "R_next"=>expm1(ll2),
        "log_C"=>logc,
        "C"=>isnothing(logc) ? nothing : exp(logc),
        "G"=>g,
        "G_x"=>lw2-2ls2,
        "G_y"=>ll2-2lw2,
        "eta_over_R_squared"=>eta_over_r2,
        "finite_step_contraction"=>ll2<ls,
        "seed_R_le_001"=>0<ls<=log1p(0.01),
        "certified"=>false,
    )
end
