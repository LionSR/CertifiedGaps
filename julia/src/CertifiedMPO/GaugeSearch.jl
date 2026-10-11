"""
    ThermalMPOGap.CertifiedMPO.GaugeSearch

Numerical gauge proposals only: Schmidt-balanced full-matrix gauges and a log-diagonal gauge
descent. Every proposal still needs a rational witness from `Norms`. Moved verbatim from
spin2 `Spin2MPOGaugeSearch`; `MN` is `CertifiedMPO.Norms`.
"""
module GaugeSearch

using LinearAlgebra
import ..Norms
const MN = Norms
export gauge_objective, optimize_gauge

import ..ChargeBlocks

"Numerical full-matrix gauges balancing left Schmidt scales against right-canonical tails."
function schmidt_gauge_proposal(cores; floor_ratio = 1e-18, physical_charges = nothing)
    if !isnothing(physical_charges)
        return ChargeBlocks.schmidt_gauges(
            cores,
            physical_charges;
            floor_ratio,
            svd_function = MN.numerical_svd,
        )
    end
    @assert 0 < floor_ratio < 1
    gauges = [ones(1, 1)]
    prefix = ones(1, 1)
    for core in cores[1:(end-1)]
        l, d, di, r = size(core)
        contracted = prefix*reshape(core, l, d*di*r)
        columns = reshape(contracted, size(prefix, 1)*d*di, r)
        prefix = Matrix(qr(columns).R)
        # Near a boundary, the auxiliary rank can exceed the physical prefix
        # dimension. Retain and regularize the full right nullspace as well.
        decomposition = MN.numerical_svd(prefix; full = true)
        singular = vcat(decomposition.S, zeros(r-length(decomposition.S)))
        scales = sqrt.(max.(singular, floor_ratio*first(decomposition.S)))
        push!(gauges, Matrix(Diagonal(scales)*decomposition.Vt))
    end
    push!(gauges, ones(1, 1))
    return gauges
end

"Numerical diagnostic only; rational local norm witnesses are still required."
function matrix_gauge_objective(cores, gauges)
    @assert length(gauges) == length(cores)+1
    value = 0.0
    for (j, core) in enumerate(cores)
        l, d, di, r = size(core)
        left = gauges[j]*reshape(core, l, d*di*r)
        transformed = reshape(reshape(left, l*d*di, r)/gauges[j+1], l, d, di, r)
        value += log(opnorm(MN.local_operator(transformed)))
    end
    return value
end

function gauge_objective(cores, x)
    dimensions = [size(c, 4) for c in cores[1:(end-1)]]
    ends = cumsum(vcat(0, dimensions))
    weights =
        vcat([zeros(1)], [x[(ends[j]+1):ends[j+1]] for j in eachindex(dimensions)], [zeros(1)])
    gradient = [zeros(length(w)) for w in weights]
    value = 0.0
    for (j, core) in enumerate(cores)
        l, d, _, r = size(core)
        transformed =
            core .* reshape(exp.(weights[j]), l, 1, 1, 1) .*
            reshape(exp.(-weights[j+1]), 1, 1, 1, r)
        W = reshape(permutedims(transformed, (2, 4, 3, 1)), d*r, d*l)
        s = MN.numerical_svd(W)
        sigma = first(s.S)
        value += log(sigma)
        cluster = findall(v->v >= sigma*(1-1e-10), s.S)
        for k in cluster
            gradient[j] += vec(sum(reshape(abs2.(s.V[:, k]), d, l); dims = 1))/length(cluster)
            gradient[j+1] -= vec(sum(reshape(abs2.(s.U[:, k]), d, r); dims = 1))/length(cluster)
        end
    end
    return value, vcat(gradient[2:(end-1)]...)
end

function optimize_gauge(cores; iterations = 200)
    x = zeros(sum(size(c, 4) for c in cores[1:(end-1)]))
    value, gradient = gauge_objective(cores, x)
    initial = value
    for iteration in 1:iterations
        norm(gradient, Inf) < 1e-7 && break
        direction = -gradient
        step = 10.0
        candidate = x+step*direction
        trial, gn = gauge_objective(cores, candidate)
        while trial > value-1e-4*step*dot(gradient, gradient) && step > 1e-8
            step /= 2
            candidate = x+step*direction
            trial, gn = gauge_objective(cores, candidate)
        end
        step <= 1e-8 && break
        x, value, gradient = candidate, trial, gn
    end
    return (; initial, final = value, x, gradient_norm = norm(gradient, Inf))
end

end
