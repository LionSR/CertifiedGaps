# Seed-scale estimate for the spin-2 Heisenberg ring, and the exact admissible region of the
# doubling recurrence. NOTHING in the scale map is a certificate: excitation sums come from a
# magnon-gas model with the numerical gap and correlation length of Todo–Kato, and the
# ring-closure energies δ_n = E₀(2n) − 2E₀(n) are differences of certified *trial* energies
# (upper bounds on E₀, so their differences are estimates). The recurrence calls are rigorous.
#
#   julia --startup-file=no --project=models/spin2 models/spin2/experiments/seed_scale.jl [out.json]
#
# Also runs with --project=julia (only the core package is used).
module SeedScale

using CertifiedGaps.Recurrence: certified_gap_bound
using CertifiedGaps.Results: write_result, run_metadata
using JSON3, SHA

export magnon_sum,
    log_purities,
    ring_closure,
    delta_model,
    seed_window,
    max_defect,
    squaring_budget,
    squaring_schedule,
    trial_energies

# Numerical reference values (Todo–Kato, PRL 87, 047203; not certified).
const GAP = 0.08917
const XI = 49.49
const VELOCITY = GAP * XI            # v = Δξ, the single-mode (relativistic) estimate

const TRIAL = joinpath(@__DIR__, "..", "results", "trial-295-symmetric-certificate.json")

"""
    trial_energies(path = TRIAL) -> Dict(N => (lower, upper))

Certified strict intervals for the D = 295 trial energy per bond, as `Float64` (display use).
"""
function trial_energies(path = TRIAL)
    data = JSON3.read(read(path, String), Dict{String,Any})
    @assert data["complete"] && data["status"] == "finite-energy-certificates-passed"
    out = Dict{Int,NamedTuple{(:lower, :upper),Tuple{Float64,Float64}}}()
    for c in data["certificates"]
        i = c["energy_per_bond_interval"]
        den = parse(BigInt, i["denominator"])
        out[c["length"]] = (
            lower = Float64(parse(BigInt, i["lower_numerator"]) // den),
            upper = Float64(parse(BigInt, i["upper_numerator"]) // den),
        )
    end
    return out
end

"""
    ring_closure(energies, n) -> δ_n

`E_trial(2n) − 2E_trial(n)` from energies per bond: an estimate of `E₀(2n) − 2E₀(n)`, not a
bound.
"""
ring_closure(energies, n) = 2n * energies[2n].upper - 2 * n * energies[n].upper

"""
    delta_model(energies; xi = nothing) -> (δ::Function, xi, amplitude, anchor)

Exponential extrapolation `δ_n = A e^{−n/ξ}`. With `xi = nothing` the decay length is fitted
through δ₁₂₈ and δ₂₅₆; otherwise the given ξ is anchored at δ₁₂₈ (a larger ξ gives the
larger, more conservative δ_n beyond 128).
"""
function delta_model(energies; xi = nothing)
    d128, d256 = ring_closure(energies, 128), ring_closure(energies, 256)
    ξ = xi === nothing ? 128 / log(d128 / d256) : float(xi)
    A = d128 * exp(128 / ξ)
    return (δ = n -> A * exp(-n / ξ), xi = ξ, amplitude = A, anchor = 128)
end

"Magnon dispersion ε(k) = √(Δ² + v²(k − π)²), minimal image in k."
dispersion(k) = sqrt(GAP^2 + VELOCITY^2 * (mod(k - π + π, 2π) - π)^2)

"""
    magnon_sum(L, β) -> Σ_{modes} 3 e^{−βε(k)}

Ideal gas of triplet magnons on the ring of `L` sites: `log(1 + R_L(β)) ≈ magnon_sum(L, β)`.
Two-magnon states enter through the exponential; interactions and bound states are ignored.
"""
magnon_sum(L, β) = 3 * sum(exp(-β * dispersion(2π * j / L)) for j in 0:(L-1))

"""
    log_purities(n, β, δn) -> (log Pₛ(n, β), log P_τ(2n, β))

From `Z_L(β) = e^{−βE₀(L)}(1 + R_L(β))` with `log(1 + R) = magnon_sum`.
"""
function log_purities(n, β, δn)
    ls = -β * δn + magnon_sum(2n, β) - 2 * magnon_sum(n, β)
    lt = magnon_sum(2n, 2β) - 2 * magnon_sum(2n, β)
    return min(ls, 0.0), min(lt, 0.0)
end

defects(n, β, δn) = map(x -> -expm1(x), log_purities(n, β, δn))

"""
    seed_window(n, βs, δn) -> (βmin, βmax, best_bound, best_β, rows)

Scan `βs`; a point is feasible when `certified_gap_bound` (rigorous for the given defects)
accepts the model defects.
"""
function seed_window(n, βs, δn)
    rows = map(βs) do β
        εs, εt = defects(n, β, δn)
        bound = (0 <= εs < 1 && 0 <= εt < 1) ? certified_gap_bound(εs, εt, β).bound : nothing
        (; β, εs, εt, bound)
    end
    feasible = filter(r -> r.bound !== nothing, rows)
    isempty(feasible) &&
        return (βmin = nothing, βmax = nothing, best_bound = nothing, best_β = nothing, rows)
    best = argmax(r -> r.bound, feasible)
    return (
        βmin = minimum(r.β for r in feasible),
        βmax = maximum(r.β for r in feasible),
        best_bound = best.bound,
        best_β = best.β,
        rows,
    )
end

admissible(εs, εt) = certified_gap_bound(εs, εt, 1).bound !== nothing

"""
    max_defect(fixed; vary = :spatial, tol = 1 // 10^6) -> (accepted, rejected)

Largest admissible defect of one kind when the other is `fixed`, by bisection on exact
rationals: `accepted` passes `certified_gap_bound`, `rejected = accepted + tol` fails.
Admissibility does not depend on β₀.
"""
function max_defect(fixed; vary = :spatial, tol = 1 // 10^6)
    test(x) = vary === :spatial ? admissible(x, fixed) : admissible(fixed, x)
    lo, hi = 0 // 1, 1 // 2
    @assert test(lo) && !test(hi)
    while hi - lo > tol
        mid = (lo + hi) // 2
        test(mid) ? (lo = mid) : (hi = mid)
    end
    return (accepted = lo, rejected = hi)
end

"""
    squaring_budget(K, final_relative_error)

If the normalised relative error at least doubles per certified squaring, the seed error plus
the sum of all compression residuals must stay below `final_relative_error / 2^K`.
"""
squaring_budget(K, final_relative_error) = final_relative_error / 2.0^K

"""
    squaring_schedule(K, β, final_relative_error) -> rows

Step-resolved budget of the squaring route: step `j` produces `β_j = β 2^{j−K}` and its
compression residual is amplified by `2^{K−j}` (forward error) or `4^{K−j}` (product of forward
and backward errors in the signed trace correction). Each row gives the residual allowed at
that step if the final error is shared equally over the `K` steps.
"""
function squaring_schedule(K, β, final_relative_error)
    [
        Dict(
            "step" => j,
            "beta" => β * 2.0^(j - K),
            "forward_budget" => final_relative_error / (K * 2.0^(K - j)),
            "signed_budget" => final_relative_error / (K * 4.0^(K - j)),
        ) for j in 1:K
    ]
end

function main(output = joinpath(@__DIR__, "seed-scale.json"))
    energies = trial_energies()
    models = Dict(
        "fit_128_256" => delta_model(energies),
        "xi_49.49_anchored_128" => delta_model(energies; xi = XI),
    )
    ns = (96, 112, 128, 144, 160, 176, 192, 208, 224, 240, 256, 288, 320, 384, 512)
    βs = 8.0:2.0:400.0
    scale = Dict{String,Any}()
    for (name, model) in models
        table = Any[]
        for n in ns
            w = seed_window(n, βs, model.δ(n))
            push!(
                table,
                Dict(
                    "n" => n,
                    "delta_n_estimate" => model.δ(n),
                    "beta_min" => w.βmin,
                    "beta_max" => w.βmax,
                    "best_gap_bound" => w.best_bound,
                    "best_beta" => w.best_β,
                ),
            )
            win =
                w.βmin === nothing ? "empty" :
                "β ∈ [$(w.βmin), $(w.βmax)], best $(w.best_bound) at β=$(w.best_β)"
            println("$name n=$n δ=$(model.δ(n)): $win")
        end
        scale[name] = Dict(
            "xi" => model.xi,
            "amplitude" => model.amplitude,
            "anchor" => model.anchor,
            "windows" => table,
        )
    end
    # Exact-diagonalisation sizes: README DMRG values (numerical), E₀(4) = −20 exact.
    small = Dict(
        "E0" => Dict(4 => -20.0, 8 => -38.52186914, 16 => -76.37583432),
        "gap" => Dict(8 => 0.54556928, 16 => 0.3109),
    )
    small_windows = map((4, 8)) do n
        δ = small["E0"][2n] - 2 * small["E0"][n]
        Dict(
            "n" => n,
            "delta_n" => δ,
            "beta_max_spatial" => 0.22 / δ,
            "beta_min_thermal" => 6 / small["gap"][2n],
            "window_nonempty" => 0.22 / δ >= 6 / small["gap"][2n],
        )
    end
    region = Dict(
        "max_spatial_defect_given_thermal" => [
            let r = max_defect(εt; vary = :spatial)
                Dict(
                    "thermal_defect" => string(εt),
                    "accepted" => string(r.accepted),
                    "rejected" => string(r.rejected),
                    "accepted_float" => Float64(r.accepted),
                )
            end for εt in (
                1 // 10^8,
                1 // 10^4,
                1 // 1000,
                1 // 100,
                3 // 100,
                5 // 100,
                8 // 100,
                1 // 10,
                12 // 100,
            )
        ],
        "max_thermal_defect_given_spatial" => [
            let r = max_defect(εs; vary = :thermal)
                Dict(
                    "spatial_defect" => string(εs),
                    "accepted" => string(r.accepted),
                    "rejected" => string(r.rejected),
                    "accepted_float" => Float64(r.accepted),
                )
            end for εs in (
                1 // 10^8,
                1 // 10^4,
                1 // 1000,
                1 // 100,
                3 // 100,
                5 // 100,
                8 // 100,
                1 // 10,
                12 // 100,
            )
        ],
    )
    for row in region["max_spatial_defect_given_thermal"]
        println("εₜ=$(row["thermal_defect"]): max εₛ ≈ $(row["accepted_float"])")
    end
    for row in region["max_thermal_defect_given_spatial"]
        println("εₛ=$(row["spatial_defect"]): max εₜ ≈ $(row["accepted_float"])")
    end
    budget = Dict(
        "squarings" => 28,
        "final_relative_error" => 0.02,
        "seed_plus_residual_budget" => squaring_budget(28, 0.02),
        "schedule_K28_degree2_seed" => squaring_schedule(28, 64.0, 0.02),
        "schedule_K16_degree12_seed" => squaring_schedule(16, 64.0, 0.02),
        "note" =>
            "relative HS error doubles per squaring; 0.02 is the slack between the " *
            "trial-based estimate of Ẑ_256(64) and the target 1.08",
    )
    result = Dict(
        "status" => "float-seed-scale-estimate",
        "certified" => false,
        "model" => Dict(
            "gap" => GAP,
            "xi" => XI,
            "velocity" => VELOCITY,
            "excitations" =>
                "ideal triplet magnon gas, ε(k)=sqrt(Δ²+v²(k−π)²), " * "log(1+R)=Σ 3e^{−βε}",
            "ring_closure" =>
                "δ_n = E_trial(2n) − 2E_trial(n) from the D=295 certificate " * "(estimate)",
        ),
        "trial_energies_per_bond_upper" => Dict(string(k) => v.upper for (k, v) in energies),
        "ring_closure_from_trial" =>
            Dict(string(n) => ring_closure(energies, n) for n in (128, 256, 512)),
        "seed_scale" => scale,
        "exact_diagonalisation_sizes" => small_windows,
        "admissible_region_exact" => region,
        "squaring_route_budget" => budget,
        "run_metadata" => run_metadata(
            Dict("ns" => collect(ns), "beta_grid" => [first(βs), step(βs), last(βs)]);
            source_sha256 = bytes2hex(sha256(read(@__FILE__))),
            status = "float estimate; admissible-region rows and recurrence calls are rigorous",
            certified = false,
        ),
    )
    write_result(output, result)
    return result
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    SeedScale.main(isempty(ARGS) ? joinpath(@__DIR__, "seed-scale.json") : only(ARGS))
end
