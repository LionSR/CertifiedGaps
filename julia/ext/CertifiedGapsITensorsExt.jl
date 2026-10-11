module ThermalMPOGapITensorsExt

using ITensors, ITensorMPS
import ThermalMPOGap.Purification: purified_log_partition

# The shared algorithm of spin2 `thermal_curve`, aklt `thermal_curve` and spin1 `thermal_logz`.
function purified_log_partition(sites, pairs, halfgate, betas; dbeta, maxdim, cutoff,
    apply_kwargs = (;), extra = nothing, progress = false, label = "", callback = nothing)
    targets = sort(unique(Float64.(betas)))
    !isempty(targets) && all(isfinite, targets) && first(targets) >= 0 ||
        throw(ArgumentError("betas must be finite and nonnegative"))
    dbeta > 0 && maxdim >= 1 && cutoff >= 0 || throw(ArgumentError("invalid evolution controls"))
    steps = round.(Int, targets ./ dbeta)
    all(isapprox.(steps .* dbeta, targets; atol = 1e-10, rtol = 1e-10)) ||
        throw(ArgumentError("each beta must be an integer multiple of dbeta"))
    n, d = length(sites), dim(first(sites))
    amplitude = MPO(sites, "Id")
    for j in 1:n
        amplitude[j] /= sqrt(d)
    end
    halfgates = [halfgate(sites, i, j, dbeta / 4) for (i, j) in pairs]
    gates = vcat(halfgates, reverse(halfgates))
    lognorm = 0.0
    result = Dict{Float64,Dict{String,Any}}()
    start = time()
    if 0 in steps
        result[0.0] = Dict("beta" => 0.0, "logZ" => n * log(float(d)), "maxlinkdim" => 1,
            "seconds" => 0.0)
    end
    for step in 1:maximum(steps)
        amplitude = apply(gates, amplitude; cutoff, maxdim, apply_kwargs...)
        growth = norm(amplitude)
        isfinite(growth) && growth > 0 || error("invalid purification norm")
        lognorm += log(growth)
        amplitude /= growth
        step in steps || continue
        beta = targets[findfirst(==(step), steps)]
        logz = n * log(float(d)) + 2lognorm
        row = Dict{String,Any}("beta" => beta, "logZ" => logz,
            "maxlinkdim" => maxlinkdim(amplitude), "seconds" => time() - start)
        extra === nothing || merge!(row, extra(logz))
        result[beta] = row
        if progress
            println("thermal $label n=$n beta=$beta dt=$dbeta chi=$maxdim logZ=$logz ",
                "seconds=$(time() - start)")
            flush(stdout)
        end
        callback === nothing || callback(result)
    end
    return result
end

end
