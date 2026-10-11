# Numerical search only: can local operator norms provide a useful MPO witness?
include(joinpath(@__DIR__, "mpo_norm_probe.jl"))

import ThermalMPOGap
using ThermalMPOGap.Results: write_result
using ThermalMPOGap.CertifiedMPO.GaugeSearch

function gauge_main(output = joinpath(@__DIR__, "../results/mpo-gauge-probe.json"))
    result = Dict{String,Any}(
        "status"=>"numerical-MPO-gauge-search",
        "certified"=>false,
        "complete"=>false,
        "uniform_gap_certified"=>false,
        "runs"=>Any[],
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    )
    n = 4
    sites = siteinds("S=2", n; conserve_qns = false)
    amplitude = MPO(sites, "Id")
    amplitude /= norm(amplitude)
    half = [exp(-0.125/4*Spin2Gap.bond_operator(sites, j, k)) for (j, k) in bonds(n)]
    gates = vcat(half, reverse(half))
    for step in 1:8
        amplitude = apply(gates, amplitude; maxdim = 64, cutoff = 1e-12)
        amplitude /= norm(amplitude)
        step in (1, 8) || continue
        cores = extract_cores(orthogonalize!(deepcopy(amplitude), 2), sites)
        started = time()
        best = optimize_gauge(cores)
        reference = opnorm(dense_operator(cores))
        row = Dict(
            "length"=>n,
            "beta"=>step/8,
            "local_norm_product_initial"=>exp(best.initial),
            "local_norm_product_optimized"=>exp(best.final),
            "independent_dense_candidate_norm"=>reference,
            "gradient_norm"=>best.gradient_norm,
            "log_diagonal_gauges"=>best.x,
            "seconds"=>time()-started,
        )
        push!(result["runs"], row)
        write_result(output, result)
        println(
            "gauge beta=$(step/8) initial=$(exp(best.initial)) optimized=$(exp(best.final)) actual=$reference seconds=$(time()-started)",
        )
        flush(stdout)
    end
    result["complete"] = true
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    gauge_main(isempty(ARGS) ? joinpath(@__DIR__, "../results/mpo-gauge-probe.json") : only(ARGS))
end
