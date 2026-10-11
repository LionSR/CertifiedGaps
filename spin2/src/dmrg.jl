mutable struct SweepHistory <: AbstractObserver
    history::Vector{Dict{String,Any}}
    maxerr::Float64
end
SweepHistory() = SweepHistory(Dict{String,Any}[], 0.0)

function ITensorMPS.measure!(obs::SweepHistory; kwargs...)
    obs.maxerr = max(obs.maxerr, kwargs[:spec].truncerr)
    if kwargs[:bond] == 1 && kwargs[:half_sweep] == 2
        push!(
            obs.history,
            Dict(
                "sweep" => kwargs[:sweep],
                "objective_energy" => real(kwargs[:energy]),
                "max_discarded_weight" => obs.maxerr,
                "maxlinkdim" => maxlinkdim(kwargs[:psi]),
            ),
        )
        obs.maxerr = 0.0
    end
    return nothing
end

function initial_state(sites, magnetization)
    n = length(sites)
    iseven(n) || throw(ArgumentError("initial state requires even n"))
    magnetization in (0, 1) || throw(ArgumentError("implemented sectors are Sz=0 and 1"))
    labels = [isodd(j) ? 1 : 5 for j in 1:n]
    magnetization == 1 && (labels[end] = 4)
    return random_mps(sites, labels; linkdims = 4)
end

function level_diagnostics(H, psi, objective, observer)
    energy = real(inner(psi', H, psi))
    variance = real(inner(H, psi, H, psi))-energy^2
    return Dict{String,Any}(
        "energy" => energy,
        "objective_energy" => real(objective),
        "variance" => variance,
        "residual_norm" => sqrt(max(0.0, variance)),
        "variance_roundoff_floor" => 100eps(Float64)*max(1, energy^2),
        "maxlinkdim" => maxlinkdim(psi),
        "sweeps" => observer.history,
        "energy_change_last_sweep" =>
            length(observer.history) < 2 ? nothing :
            abs(
                observer.history[end]["objective_energy"]-observer.history[end-1]["objective_energy"],
            ),
    )
end

"""Finite periodic spectrum estimates, including the first Sz=0 excitation.

The Sz=1 sector energy alone can miss a low singlet. Both channels are
reported, with residuals and overlaps; their difference is NOT a rigorous
gap lower bound and does not certify the two-purity seeds.
"""
function dmrg_levels(
    n;
    maxdim = 128,
    nsweeps = 10,
    cutoff = 1e-10,
    seed = 20261007,
    excited = true,
    weight = 10.0,
    outputlevel = 1,
    checkpoint_dir = nothing,
)
    bonds(n) # validate periodic size
    maxdim >= 4 || throw(ArgumentError("maxdim must be at least four"))
    nsweeps >= 4 || throw(ArgumentError("use at least four sweeps"))
    Random.seed!(seed)
    sites = siteinds("S=2", n; conserve_qns = true)
    H = hamiltonian(sites)
    schedule = [min(maxdim, d) for d in (16, 32, 64, 128, 256, 512)]
    # Reach the requested dimension by the middle of the run.
    ramp = [schedule[min(j, length(schedule))] for j in 1:div(nsweeps, 2)]
    dims = vcat(ramp, fill(maxdim, nsweeps-length(ramp)))
    noise = vcat([1e-6, 1e-7], fill(0.0, nsweeps-2))
    levels = Dict{String,Any}()
    ground = nothing
    for (label, m, penalize) in
        [("ground_sz0", 0, false), ("lowest_sz1", 1, false), ("excited_sz0", 0, true)]
        penalize && !excited && continue
        observer = SweepHistory()
        psi = initial_state(sites, m)
        start = time()
        if penalize
            objective, psi = dmrg(
                H,
                [ground],
                psi;
                nsweeps,
                maxdim = dims,
                cutoff,
                noise,
                observer,
                weight,
                outputlevel,
                eigsolve_krylovdim = 8,
                eigsolve_tol = 1e-11,
            )
        else
            objective, psi = dmrg(
                H,
                psi;
                nsweeps,
                maxdim = dims,
                cutoff,
                noise,
                observer,
                outputlevel,
                eigsolve_krylovdim = 8,
                eigsolve_tol = 1e-11,
            )
        end
        normalize!(psi)
        levels[label] = level_diagnostics(H, psi, objective, observer)
        levels[label]["seconds"] = time()-start
        penalize && (levels[label]["ground_overlap"] = abs(inner(ground, psi)))
        if !isnothing(checkpoint_dir)
            mkpath(checkpoint_dir)
            path = joinpath(checkpoint_dir, "n$(n)-chi$(maxdim)-seed$(seed)-$(label).jls")
            serialize(path, (; psi, n, maxdim, seed, label, energy = levels[label]["energy"]))
            levels[label]["checkpoint"] = abspath(path)
            levels[label]["checkpoint_sha256"] = bytes2hex(sha256(read(path)))
        end
        label == "ground_sz0" && (ground = psi)
        println(
            "n=$n chi=$maxdim $label E=$(levels[label]["energy"]) variance=$(levels[label]["variance"])",
        )
        flush(stdout)
    end
    ground_energy = levels["ground_sz0"]["energy"]
    return Dict(
        "n" => n,
        "boundary" => "periodic",
        "maxdim" => maxdim,
        "seed" => seed,
        "levels" => levels,
        "sz1_gap_estimate" => levels["lowest_sz1"]["energy"]-ground_energy,
        "sz0_gap_estimate" => excited ? levels["excited_sz0"]["energy"]-ground_energy : nothing,
        "certified" => false,
    )
end
