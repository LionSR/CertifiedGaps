"""
    CertifiedGaps.CertifiedMPO.PositiveBounds

Exact entrywise majorants M ≤ U V with nonnegative rational factors, from numerical NMF and
covering-LP candidates; acceptance is an exact matrix inequality. Moved verbatim from spin2
`Spin2PositiveBounds`.
"""
module PositiveBounds

using LinearAlgebra
import Nemo
using ..CertifiedMPO: F

using Random: MersenneTwister

"Numerical covering LP through its feasible dual; exact acceptance happens later."
function covering_candidate(V, b; weights = b, tolerance = 1e-11, maxpivots = 3000)
    rank, columns = size(V)
    @assert length(b) == columns && length(weights) == columns
    @assert all(>=(0), V) && all(>=(0), b) && all(>=(0), weights)
    scale = maximum(b)
    scale == 0 && return zeros(rank)
    cost = V*(weights/maximum(weights))
    maximum(cost) > 0 || return zeros(rank)
    active = findall(>(0), cost)
    active_rank = length(active)
    # Dual: max b'y with V*y <= cost, y >= 0. Its slack basis is feasible.
    # Equilibrate every constraint to avoid very small initial slack values.
    tableau =
        hcat(Float64.(V[active, :]) ./ cost[active], Matrix{Float64}(I, active_rank, active_rank))
    rhs = ones(active_rank)
    reduced = vcat(Float64.(b)/scale, zeros(active_rank))
    basis = collect((columns+1):(columns+active_rank))
    for _ in 1:maxpivots
        entering = findfirst(>(tolerance), reduced)
        isnothing(entering) && break
        eligible = findall(>(tolerance), tableau[:, entering])
        isempty(eligible) && break
        ratios = rhs[eligible] ./ tableau[eligible, entering]
        minimum_ratio = minimum(ratios)
        ties = [
            eligible[i] for i in eachindex(eligible) if
            ratios[i] <= minimum_ratio+tolerance*max(1, abs(minimum_ratio))
        ]
        leaving = ties[argmin(basis[ties])]
        pivot = tableau[leaving, entering]
        tableau[leaving, :] /= pivot
        rhs[leaving] /= pivot
        for i in 1:active_rank
            i == leaving && continue
            multiplier = tableau[i, entering]
            tableau[i, :] -= multiplier*tableau[leaving, :]
            rhs[i] -= multiplier*rhs[leaving]
        end
        multiplier = reduced[entering]
        reduced -= multiplier*tableau[leaving, :]
        basis[leaving] = entering
    end
    result = zeros(rank)
    result[active] = scale*max.(-reduced[(columns+1):end], 0.0) ./ cost[active]
    return result
end

function improve_cover(A, U, V; sweeps = 2)
    for _ in 1:sweeps
        for i in axes(A, 1)
            U[i, :] = covering_candidate(V, A[i, :])
        end
        for j in axes(A, 2)
            V[:, j] = covering_candidate(transpose(U), A[:, j])
        end
        for k in axes(U, 2)
            a, b = norm(U[:, k]), norm(V[k, :])
            if a > 0 && b > 0
                f = sqrt(b/a)
                U[:, k] *= f
                V[k, :] /= f
            end
        end
    end
    return U, V
end

"Numerical nonnegative factorization; its output is only a candidate."
function nmf_candidate(M, rank; iterations = 1000, seed = 17, underweight = 1.0)
    all(>=(0), M) || throw(ArgumentError("matrix must be nonnegative"))
    1 <= rank <= min(size(M)...) || throw(ArgumentError("invalid rank"))
    scale = maximum(M)
    scale == 0 && return zeros(size(M, 1), rank), zeros(rank, size(M, 2))
    A = Float64.(M)/Float64(scale)
    s = svd(A)
    U, V = zeros(size(A, 1), rank), zeros(rank, size(A, 2))
    U[:, 1] = sqrt(s.S[1])*abs.(s.U[:, 1])
    V[1, :] = sqrt(s.S[1])*abs.(s.Vt[1, :])
    for k in 2:rank
        a, b = s.U[:, k], s.Vt[k, :]
        ap, am, bp, bm = max.(a, 0), max.(-a, 0), max.(b, 0), max.(-b, 0)
        positive, negative = norm(ap)*norm(bp), norm(am)*norm(bm)
        x, y, mass = positive >= negative ? (ap, bp, positive) : (am, bm, negative)
        if mass > 0
            U[:, k] = sqrt(s.S[k]*mass)*x/norm(x)
            V[k, :] = sqrt(s.S[k]*mass)*y/norm(y)
        end
    end
    rng = MersenneTwister(seed)
    U += 1e-6*rand(rng, size(U)...)
    V += 1e-6*rand(rng, size(V)...)
    for iteration in 1:iterations
        if underweight == 1
            U .*= (A*transpose(V)) ./ max.(U*(V*transpose(V)), 1e-300)
            V .*= (transpose(U)*A) ./ max.((transpose(U)*U)*V, 1e-300)
        else
            penalty = underweight^(iteration/iterations)
            approximation = U*V
            weights = ifelse.(approximation .< A, penalty, 1.0)
            U .*=
                ((weights .* A)*transpose(V)) ./
                max.((weights .* approximation)*transpose(V), 1e-300)
            approximation = U*V
            weights = ifelse.(approximation .< A, penalty, 1.0)
            V .*=
                (transpose(U)*(weights .* A)) ./
                max.(transpose(U)*(weights .* approximation), 1e-300)
        end
        if iteration % 50 == 0
            for k in 1:rank
                a, b = norm(U[:, k]), norm(V[k, :])
                if a > 0 && b > 0
                    f = sqrt(b/a)
                    U[:, k] *= f
                    V[k, :] /= f
                end
            end
        end
    end
    return sqrt(Float64(scale))*U, sqrt(Float64(scale))*V
end

"Rational bracket whose validity is checked by outward-rounded comparison."
function rational_bracket(x; bits = 80)
    denominator = big(2)^bits
    center = round(BigInt, BigFloat(x)*BigFloat(denominator))
    radius = big(2)
    RR = parent(x)
    while !(x > RR((center-radius)//denominator) && x < RR((center+radius)//denominator))
        radius *= 2
    end
    return (center-radius)//denominator, (center+radius)//denominator
end

"Split off troublesome residual columns; all selected envelopes are exact."
function split_envelope(residual, base, left, right, rank; bits = 40)
    A = Float64.(residual) .* left .* transpose(right)
    B = Float64.(base) .* left .* transpose(right)
    groups = [collect(axes(A, 2))]
    function score(cols)
        upper = maximum(A[:, cols]; dims = 2)
        return sum(abs2, B[:, cols] .+ upper)
    end
    for _ in 2:min(rank, size(A, 2))
        gain, selected, column = -Inf, 0, 0
        for (i, cols) in enumerate(groups)
            length(cols) > 1 || continue
            before = score(cols)
            for j in cols
                after = score([j]) + score(filter(!=(j), cols))
                if before-after > gain
                    gain, selected, column = before-after, i, j
                end
            end
        end
        selected == 0 && break
        filter!(!=(column), groups[selected])
        push!(groups, [column])
    end
    denominator = big(2)^bits
    v = [max(big(1), round(BigInt, denominator/r))//denominator for r in right]
    U, V = zeros(F, size(A, 1), length(groups)), zeros(F, length(groups), size(A, 2))
    for (k, cols) in enumerate(groups)
        V[k, cols] = v[cols]
        for i in axes(A, 1)
            U[i, k] = maximum(residual[i, j]/v[j] for j in cols)
        end
    end
    @assert all(U*V .>= residual)
    return U, V
end

"""
Exact entrywise majorant M <= U V with nonnegative rational factors.

Quantize a numerical rank-r factorization and add a rank-one upper envelope
of its positive residual. Acceptance is an exact rational matrix inequality,
independent of convergence or correctness of the numerical optimizer.
"""
function upper_factorization(
    M::AbstractMatrix{F},
    rank;
    iterations = 1000,
    bits = 40,
    row_metric = ones(size(M, 1)),
    col_metric = ones(size(M, 2)),
    underweight = 1.0,
    cover_sweeps = 0,
    correction_rank = 1,
)
    all(>=(0), M) || throw(ArgumentError("matrix must be nonnegative"))
    1 <= rank <= min(size(M)...) || throw(ArgumentError("invalid rank"))
    if rank == min(size(M)...)
        U, V =
            size(M, 1) <= size(M, 2) ? (Matrix{F}(I, size(M, 1), size(M, 1)), copy(M)) :
            (copy(M), Matrix{F}(I, size(M, 2), size(M, 2)))
        return (; U, V, bound = copy(M), rank, excess = zero(F), mass = sum(M))
    end
    left, right = sqrt.(Float64.(row_metric)), sqrt.(Float64.(col_metric))
    @assert all(>(0), left) && all(>(0), right)
    U0, V0 = nmf_candidate(Float64.(M) .* left .* transpose(right), rank; iterations, underweight)
    baseline = enclose_candidate(M, U0, V0, left, right; bits, correction_rank)
    if cover_sweeps > 0
        U1, V1 = improve_cover(
            Float64.(M) .* left .* transpose(right),
            copy(U0),
            copy(V0);
            sweeps = cover_sweeps,
        )
        if all(isfinite, U1) && all(isfinite, V1)
            covered = enclose_candidate(M, U1, V1, left, right; bits, correction_rank)
            score(c) = sum(abs2, Float64.(c.bound) .* left .* transpose(right))
            score(covered) < score(baseline) && return covered
        end
    end
    return baseline
end

function enclose_candidate(M, U0, V0, left, right; bits, correction_rank)
    @assert correction_rank >= 1
    U0, V0 = copy(U0), copy(V0)
    U0 ./= left
    V0 ./= transpose(right)
    denominator = big(2)^bits
    U = [round(BigInt, u*denominator)//denominator for u in U0]
    V = [round(BigInt, v*denominator)//denominator for v in V0]
    residual = max.(M-U*V, zero(F))
    if !iszero(residual)
        # Compare row and column envelopes. Both are exact majorants.
        ur, vr = maximum(residual; dims = 2), ones(F, 1, size(M, 2))
        uc, vc = ones(F, size(M, 1), 1), maximum(residual; dims = 1)
        weighted_v = [round(BigInt, denominator/r)//denominator for r in right]
        weighted_u = [round(BigInt, denominator/l)//denominator for l in left]
        @assert all(>(0), weighted_u) && all(>(0), weighted_v)
        u1 = [maximum(residual[i, j]/weighted_v[j] for j in axes(M, 2)) for i in axes(M, 1)]
        v1 = [maximum(residual[i, j]/weighted_u[i] for i in axes(M, 1)) for j in axes(M, 2)]
        envelopes = [
            (ur, vr),
            (uc, vc),
            (reshape(u1, :, 1), reshape(weighted_v, 1, :)),
            (reshape(weighted_u, :, 1), reshape(v1, 1, :)),
        ]
        if correction_rank > 1
            push!(envelopes, split_envelope(residual, U*V, left, right, correction_rank; bits))
            a, b = split_envelope(
                transpose(residual),
                transpose(U*V),
                right,
                left,
                correction_rank;
                bits,
            )
            push!(envelopes, (Matrix(transpose(b)), Matrix(transpose(a))))
        end
        base = Float64.(U*V)
        scores = [
            sum(abs2, (base+Float64.(a)*Float64.(b)) .* left .* transpose(right)) for
            (a, b) in envelopes
        ]
        u, v = envelopes[argmin(scores)]
        U, V = hcat(U, u), vcat(V, v)
    end
    @assert all(>=(0), U) && all(>=(0), V)
    bound = U*V
    @assert all(bound .>= M)
    return (; U, V, bound, rank = size(U, 2), excess = sum(bound-M), mass = sum(M))
end

end
