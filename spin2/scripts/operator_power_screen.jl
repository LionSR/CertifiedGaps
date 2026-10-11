"Numerically select three independent gauge powers; no screening value is a certificate."
function operator_power_screen(
    input,
    square,
    proof;
    powers = 0.0:0.05:0.6,
    precision = 192,
    second_input = input,
)
    powers = collect(powers)
    @assert !isempty(powers) && 0.5 in powers
    n, d = length(input), size(first(input), 2)
    weights = square.weights
    physical =
        isnothing(weights) ? ones(d, d) :
        sqrt.([Float64(weights[t]/weights[s]) for s in 1:d, t in 1:d])
    families = (
        (OE.first_input(input, square), proof.first_gauges),
        (second_input, proof.second_gauges),
        (square.cores, proof.prefix_gauges),
    )
    scales, loglocal = Any[], Any[]
    for (cores, gauges) in families
        rows = [
            g isa AbstractVector ? abs.(Float64.(g)) : vec(sqrt.(sum(abs2, Float64.(g); dims = 2))) for g in gauges
        ]
        @assert all(all(>(0), r) for r in rows)
        family_scales = [[r .^ (2power-1) for r in rows] for power in powers]
        @assert all(all(all(x -> isfinite(x) && x > 0, r) for r in s) for s in family_scales)
        push!(scales, family_scales)
        transformed = MN.gauge_transform(cores, gauges; interval_precision = precision)
        numeric = [Float64.(c) .* reshape(physical, 1, d, d, 1) for c in transformed]
        family_logs = [zeros(n) for _ in powers]
        for (k, s) in enumerate(family_scales), j in 1:n
            changed = numeric[j] .* reshape(s[j], :, 1, 1, 1) ./ reshape(s[j+1], 1, 1, 1, :)
            family_logs[k][j] = log(opnorm(MN.local_operator(changed)))
        end
        push!(loglocal, family_logs)
    end
    # Contract the prefix index once per prefix power. The remaining two
    # positive-vector contractions make the three-parameter search inexpensive.
    buckets, logsizes = Any[], Float64[]
    OE.foreach_residual(
        input,
        square,
        proof.first_gauges,
        proof.second_gauges,
        proof.prefix_gauges;
        arithmetic = :interval,
        precision,
        second_input,
    ) do j, e, _, _
        q, _, _, ab = size(e)
        r, rb = size(input[j], 4), size(second_input[j], 4)
        @assert ab == r*rb
        numeric = Float64.(e) .* reshape(physical, 1, d, d, 1)
        largest = maximum(abs, numeric)
        push!(logsizes, iszero(largest) ? -Inf : log(largest))
        normalized = iszero(largest) ? numeric : numeric/largest
        bucket = reshape(sum(abs2, reshape(normalized, q, d, d, r, rb); dims = (2, 3)), q, r*rb)
        push!(
            buckets,
            [reshape(transpose(scales[3][k][j] .^ 2)*bucket, r, rb) for k in eachindex(powers)],
        )
        if j % 16 == 0 || j == n-1
            println("Operator gauge screen: residual $j/$(n-1)")
            flush(stdout)
        end
    end
    prefix_logs = [vcat(0.0, cumsum(v)) for v in loglocal[3]]
    tails = [[vcat(reverse(cumsum(reverse(v))), 0.0) for v in family] for family in loglocal[1:2]]
    function objective(indices)
        a, b, c = indices
        terms = Float64[]
        for j in 1:(n-1)
            u, v = inv.(scales[1][a][j+1] .^ 2), inv.(scales[2][b][j+1] .^ 2)
            norm2 = dot(u, buckets[j][c]*v)
            @assert norm2 >= 0 && isfinite(norm2)
            term = logsizes[j]+0.5log(norm2)+prefix_logs[c][j]+tails[1][a][j+1]+tails[2][b][j+1]
            term == -Inf && continue
            isfinite(term) || return Inf
            push!(terms, term)
        end
        isempty(terms) && return -Inf
        largest = maximum(terms)
        return largest+log(sum(exp.(terms .- largest)))
    end
    baseline = findfirst(==(0.5), powers)
    best = (baseline, baseline, baseline)
    best_log = objective(best)
    baseline_log = best_log
    for a in eachindex(powers), b in eachindex(powers), c in eachindex(powers)
        value = objective((a, b, c))
        if value < best_log
            best, best_log = (a, b, c), value
        end
    end
    @assert isfinite(baseline_log) && isfinite(best_log)
    println(
        "Proposed operator-error powers: $(powers[collect(best)]); estimated improvement=$(exp(baseline_log-best_log))",
    )
    flush(stdout)
    selected = ntuple(k -> [MN.F.(s) for s in scales[k][best[k]]], 3)
    return (; powers = powers[collect(best)], scales = selected, baseline_log, best_log)
end
