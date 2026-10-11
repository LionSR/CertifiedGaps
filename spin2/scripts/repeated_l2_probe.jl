# Dense reference for certified Hilbert--Schmidt error propagation. This is not
# yet an MPO algorithm: its operator-norm witness enumerates four-site sectors.
include(joinpath(@__DIR__, "positive_compression_probe.jl"))
using ThermalMPOGap.Results: write_result

const SECTORS =
    [[code for code in 0:624 if sum(mod(div(code, 5^j), 5)-2 for j in 0:3) == mag] for mag in -8:8]
const METRIC = [local_metric(i) for i in 0:624]
location(r, c) = (mod(r, 25)+25mod(c, 25)+1, div(r, 25)+25div(c, 25)+1)
weight(code) = prod(binomial(4, mod(div(code, 5^j), 5)) for j in 0:3)
norm2(M) = sum(M[i, j]^2*METRIC[i]*METRIC[j] for i in 1:625, j in 1:625)

function hermitize(M)
    B = zero(M)
    for states in SECTORS, r in states, c in states
        i, j = location(r, c)
        k, l = location(c, r)
        B[i, j] = (M[i, j] + F(weight(r), weight(c))*M[k, l])/2
    end
    return B
end

function squared_interval(M)
    RR = Nemo.ArbField(192)
    lower, upper = zero(M), zero(M)
    for states in SECTORS
        block = [M[location(r, c)...] for r in states, c in states]
        squared = Nemo.matrix(RR, block)^2
        for (i, r) in enumerate(states), (j, c) in enumerate(states)
            lower[location(r, c)...], upper[location(r, c)...] = PB.rational_bracket(squared[i, j])
        end
    end
    return hermitize(lower), hermitize(upper)
end

function operator_norm_upper(M; details = false)
    # A is Hermitian in the physical metric. For each strictly positive v,
    # ||A|| <= rho(|A|) <= max_i (|A|v)_i/v_i. Every final comparison is rational.
    @assert hermitize(M) == M
    result = zero(F)
    vectors = Vector{F}[]
    for states in SECTORS
        block = abs.([M[location(r, c)...] for r in states, c in states])
        if iszero(block)
            push!(vectors, ones(F, length(states)))
            continue
        end
        numerical = eigen(Float64.(block))
        i = argmax(real.(numerical.values))
        v0 = abs.(real.(numerical.vectors[:, i]))
        v0 /= maximum(v0)
        denominator = big(2)^40
        v = [max(big(1), round(BigInt, (a+1e-10)*denominator))//denominator for a in v0]
        push!(vectors, v)
        result = max(result, maximum((block*v) ./ v))
    end
    return details ? (; upper = result, vectors) : result
end

function verify_operator_norm_witness(M, proof)
    @assert hermitize(M) == M && proof.upper >= 0
    @assert length(proof.vectors) == length(SECTORS)
    for (states, v) in zip(SECTORS, proof.vectors)
        @assert length(states) == length(v) && all(>(0), v)
        block=abs.([M[location(r, c)...] for r in states, c in states])
        @assert all(block*v .<= proof.upper*v)
    end
    return true
end

function rounded_svd(block, left_metric, right_metric, k, bits)
    k == min(size(block)...) && return copy(block)
    k == 0 && return zero(block)
    denominator = big(2)^bits
    left, right = sqrt.(Float64.(left_metric)), sqrt.(Float64.(right_metric))
    s = svd(Float64.(block) .* left .* transpose(right))
    u = (s.U[:, 1:k]*Diagonal(sqrt.(s.S[1:k]))) ./ left
    v = (Diagonal(sqrt.(s.S[1:k]))*s.Vt[1:k, :]) ./ transpose(right)
    U = [round(BigInt, a*denominator)//denominator for a in u]
    V = [round(BigInt, a*denominator)//denominator for a in v]
    return U*V
end

swap_local(i) = div(i-1, 25)+25mod(i-1, 25)+1

function parity_blocks(M, indices)
    lookup = Dict(i=>k for (k, i) in enumerate(indices))
    blocks = Any[]
    for sign in (1, -1)
        basis = Any[]
        for i in indices
            j = swap_local(i)
            j < i && continue
            i == j && sign == -1 && continue
            support = i == j ? [lookup[i]] : [lookup[i], lookup[j]]
            coefficients = i == j ? F[1] : F[1, sign*METRIC[i]]
            squared_norm = i == j ? F(1) : 2METRIC[i]
            push!(basis, (; support, coefficients, squared_norm))
        end
        C = zeros(F, length(basis), length(basis))
        for (a, u) in enumerate(basis), (b, v) in enumerate(basis)
            C[a, b] =
                sum(
                    c*d*METRIC[indices[i]]*METRIC[indices[j]]*M[i, j] for
                    (i, c) in zip(u.support, u.coefficients) for
                    (j, d) in zip(v.support, v.coefficients)
                )/(u.squared_norm*v.squared_norm)
        end
        metric = getproperty.(basis, :squared_norm)
        values =
            svdvals(Float64.(C) .* sqrt.(Float64.(metric)) .* transpose(sqrt.(Float64.(metric))))
        push!(blocks, (; basis, C, metric, values))
    end
    return blocks
end

function svd_candidate(M, rankcap; bits = 40)
    @assert hermitize(M) == M
    B = zero(M)
    rank = 0
    groups = Dict(q=>findall(i->local_charge(i-1)==q, 1:625) for q in -8:8)
    for q in 0:8
        rows, cols = groups[q], groups[-q]
        block = M[rows, cols]
        k = min(rankcap, size(block)...)
        if q > 0
            C = rounded_svd(block, METRIC[rows], METRIC[cols], k, bits)
            B[rows, cols] = C
            B[swap_local.(rows), swap_local.(cols)] = C .* METRIC[rows] .* transpose(METRIC[cols])
            rank += 2k
        elseif k == min(size(block)...)
            B[rows, cols] = block
            rank += k
        else
            blocks = parity_blocks(block, rows)
            singular = sort(
                [(v, p) for (p, b) in enumerate(blocks) for v in b.values];
                by = first,
                rev = true,
            )
            for (p, b) in enumerate(blocks)
                keep = count(x->last(x)==p, singular[1:k])
                C = rounded_svd(b.C, b.metric, b.metric, keep, bits)
                for (a, u) in enumerate(b.basis), (c, v) in enumerate(b.basis)
                    for (i, x) in zip(u.support, u.coefficients),
                        (j, y) in zip(v.support, v.coefficients)

                        B[rows[i], cols[j]] += x*y*C[a, c]
                    end
                end
                rank += keep
            end
        end
    end
    @assert hermitize(B) == B
    return B, rank
end

function sqrt_upper(x)
    x == 0 && return zero(F)
    return last(PB.rational_bracket(sqrt(Nemo.ArbField(192)(x))))
end

function run_l2(output; rankcaps = (16, 32, 64), beta_target = 4, doublings = 5)
    result = Dict{String,Any}(
        "status"=>"four-site-certified-hilbert-schmidt-compression-probe",
        "complete"=>false,
        "uniform_gap_certified"=>false,
        "length"=>4,
        "beta_target"=>beta_target,
        "doublings"=>doublings,
        "steps_final"=>2^(doublings+1),
        "runs"=>Any[],
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "reference_runner_sha256"=>bytes2hex(
            sha256(read(joinpath(@__DIR__, "positive_compression_probe.jl"))),
        ),
        "source_sha256"=>Dict(
            name=>bytes2hex(sha256(read(joinpath(@__DIR__, "../src", name)))) for
            name in ("ThermalBounds.jl", "../../julia/src/CertifiedMPO/PositiveBounds.jl")
        ),
    )
    references = [amplitude_interval(beta_target//2^(doublings-k), 2^(k+1)) for k in 0:doublings]
    for rankcap in rankcaps
        lower, upper, _ = first(references)
        M = hermitize((lower+upper)/2)
        error = sqrt_upper(norm2((upper-lower)/2))
        started = time()
        for k in 0:doublings
            norm_bound = operator_norm_upper(M)
            if k > 0
                lower, upper = squared_interval(M)
                M = (lower+upper)/2
                rounding = sqrt_upper(norm2((upper-lower)/2))
                error = error*(2norm_bound+error)+rounding
            end
            candidate, rank = svd_candidate(M, rankcap)
            truncation = sqrt_upper(norm2(candidate-M))
            error += truncation
            # Keep the propagated rational endpoint small in representation size.
            grid = big(2)^80
            error = cld(numerator(error)*grid, denominator(error))//grid
            M = candidate
            lo_norm, hi_norm = PB.rational_bracket(sqrt(Nemo.ArbField(192)(norm2(M))))
            trace_lower = max(lo_norm-error, zero(F))^2
            trace_upper = (hi_norm+error)^2
            reference = last(references[k+1])
            @assert reference > trace_lower && reference < trace_upper
            row = Dict(
                "rankcap_per_charge"=>rankcap,
                "beta"=>string(beta_target//2^(doublings-k)),
                "doublings_completed"=>k,
                "represented_rank_upper"=>rank,
                "certified_error_upper"=>string(error),
                "operator_norm_upper_before_step"=>string(norm_bound),
                "thermal_trace_lower"=>string(trace_lower),
                "thermal_trace_upper"=>string(trace_upper),
                "direct_trotter_trace_ball"=>string(reference),
                "cumulative_inflation_display"=>Float64(trace_upper)/Float64(reference),
                "error_upper_display"=>Float64(error),
                "seconds"=>time()-started,
            )
            push!(result["runs"], row)
            write_result(output, result)
            println(
                "HS compression cap=$rankcap stage=$k beta=$(row["beta"]) upper/reference=$(row["cumulative_inflation_display"]) error=$(Float64(error)) seconds=$(time()-started)",
            )
            flush(stdout)
        end
    end
    result["complete"] = true
    result["completed_utc"] = string(now(UTC))
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_l2(
        isempty(ARGS) ? joinpath(@__DIR__, "../results/repeated-l2-compression.json") : only(ARGS),
    )
end
