"""
    CertifiedGaps.CertifiedMPO.MPOAlgebra

Exact rational MPO algebra: direct sums, local products, exact right-tail compression with
rational intertwiners and its replay, invertible conditioning gauges, the certified initial
HS error between two MPOs, and sparse expansion for small controls.

Moved verbatim from the d-generic part of spin2 `Spin2PolynomialMPO`. Where spin2 hard-coded
the spin-2 local dimension 5 or the weights `WEIGHTS = [1, 4, 6, 4, 1]`, the functions take
`d` or `weights` keywords whose defaults reproduce spin2 exactly:

* `identity_mpo(n; d = 5)`;
* `linear_combination(mpos, coefficients; d = size(first(first(mpos)), 2))` (spin2 always
  passed d = 5 cores, and any other d failed there);
* `condition_right_exact(cores; weights = WEIGHTS, …)` and
  `initial_error(exact, candidate; weights = WEIGHTS, …)`.
"""
module MPOAlgebra

using LinearAlgebra, SparseArrays
import Nemo
using ..CertifiedMPO: nemo_array
import ..Norms, ..Squaring, ..PolynomialBounds
const MS = Squaring
const MN = Norms
const PG = PolynomialBounds
const F = MN.F

"Binomial physical weights of spin2's polynomial basis (S = 2); the default `weights`."
const WEIGHTS = F[1, 4, 6, 4, 1]

identity_mpo(n; d = 5) = fill(reshape(Matrix{F}(I, d, d), 1, d, d, 1), n)

"Exact direct sum of MPOs; coefficients enter only the first core."
function linear_combination(mpos, coefficients; d = size(first(first(mpos)), 2))
    n = length(first(mpos))
    @assert n >= 2 && length(mpos) == length(coefficients)
    @assert all(length(m)==n for m in mpos)
    result = Array{F,4}[]
    cache = Dict{Tuple,Array{F,4}}()
    for j in 1:n
        key = (j==1, j==n, Tuple(objectid(m[j]) for m in mpos))
        if haskey(cache, key)
            push!(result, cache[key])
            continue
        end
        left = j==1 ? 1 : sum(size(m[j], 1) for m in mpos)
        right = j==n ? 1 : sum(size(m[j], 4) for m in mpos)
        core = zeros(F, left, d, d, right)
        lo, ro = 0, 0
        for (m, c) in zip(mpos, coefficients)
            l, _, _, r = size(m[j])
            lr, rr = (j==1 ? (1:1) : ((lo+1):(lo+l))), (j==n ? (1:1) : ((ro+1):(ro+r)))
            core[lr, :, :, rr] += (j==1 ? c : F(1))*m[j]
            lo += l
            ro += r
        end
        cache[key] = core
        push!(result, core)
    end
    return result
end

"Exact local multiplication, exploiting sparse nonzero auxiliary blocks."
function product_core(A, B)
    la, d, mid, ra = size(A)
    lb, mid2, di, rb = size(B)
    @assert mid == mid2
    C = zeros(F, la*lb, d, di, ra*rb)
    blocks(X) =
        [(l, r, X[l, :, :, r]) for l in axes(X, 1), r in axes(X, 4) if any(!iszero, X[l, :, :, r])]
    for (l1, r1, a) in blocks(A), (l2, r2, b) in blocks(B)
        C[l1+la*(l2-1), :, :, r1+ra*(r2-1)] += a*b
    end
    return C
end

function product_mpo(A, B)
    @assert length(A) == length(B)
    cache = Dict{Tuple{UInt,UInt},Array{F,4}}()
    return [get!(cache, (objectid(a), objectid(b))) do
        product_core(a, b)
    end for (a, b) in zip(A, B)]
end

"Remove exact right-tail linear dependencies, with rational local intertwiners."
function exact_compress_right(cores)
    n = length(cores)
    @assert n >= 2
    compressed = copy(cores)
    transforms = [zeros(F, 0, 0) for _ in 1:(n+1)]
    transforms[1] = transforms[end] = ones(F, 1, 1)
    for j in n:-1:2
        l, d, di, r = size(compressed[j])
        M = Nemo.matrix(Nemo.QQ, reshape(compressed[j], l, d*di*r))
        rank, reduced = Nemo.rref(M)
        @assert rank > 0
        basis = reduced[1:rank, 1:Nemo.ncols(M)]
        pivots = [findfirst(k -> !iszero(basis[i, k]), 1:Nemo.ncols(M)) for i in 1:rank]
        L = M[1:l, pivots]
        @assert M == L*basis
        transforms[j] = nemo_array(L)
        compressed[j] = reshape(nemo_array(basis), rank, d, di, r)
        lp, dp, dip, rp = size(compressed[j-1])
        @assert rp == l
        previous = Nemo.matrix(Nemo.QQ, reshape(compressed[j-1], lp*dp*dip, rp))*L
        compressed[j-1] = reshape(nemo_array(previous), lp, dp, dip, rank)
    end
    @assert verify_exact_compression(cores, compressed, transforms)
    return (; cores = compressed, transforms)
end

"Check local intertwining identities; they telescope to exact MPO equality."
function verify_exact_compression(original, compressed, transforms)
    n = length(original)
    @assert length(compressed) == n && length(transforms) == n+1
    @assert transforms[1] == transforms[end] == ones(F, 1, 1)
    for j in 1:n
        l, d, di, r = size(original[j])
        lc, dc, dic, rc = size(compressed[j])
        @assert (d, di) == (dc, dic)
        @assert size(transforms[j]) == (l, lc) && size(transforms[j+1]) == (r, rc)
        rhs =
            Nemo.matrix(Nemo.QQ, transforms[j]) *
            Nemo.matrix(Nemo.QQ, reshape(compressed[j], lc, d*di*rc))
        lhs_flat =
            Nemo.matrix(Nemo.QQ, reshape(original[j], l*d*di, r)) *
            Nemo.matrix(Nemo.QQ, transforms[j+1])
        for a in 1:rc, t in 1:di, s in 1:d, b in 1:l
            @assert lhs_flat[b+l*(s-1)+l*d*(t-1), a] == rhs[b, s+d*(t-1)+d*di*(a-1)]
        end
    end
    return true
end

encode_transforms(transforms) =
    [Dict("dimensions"=>collect(size(T)), "entries"=>string.(vec(T))) for T in transforms]

decode_transforms(encoded) =
    [reshape(MN.parse_fraction.(T["entries"]), Int.(T["dimensions"])...) for T in encoded]

"Condition right-tail bases by invertible rational gauges, preserving exact equality."
function condition_right_exact(cores; precision = 256, bits = 80, weights = WEIGHTS)
    result = copy(cores)
    transforms = [ones(F, 1, 1) for _ in 1:(length(cores)+1)]
    for j in length(cores):-1:2
        l, d, di, r = size(result[j])
        M = reshape(result[j], l, d*di*r)
        metric = repeat(vec([F(weights[t])/F(weights[s]) for s in 1:d, t in 1:di]), r)
        gram, _ = MN.weighted_gram(Matrix(transpose(M)); row_metric = metric, precision)
        @assert size(M, 2) >= l
        L = Nemo.cholesky(gram)
        T = setprecision(BigFloat, precision) do
            proposed = [i>=k ? BigFloat(L[i, k]) : BigFloat(0) for i in 1:l, k in 1:l]
            rounded = MS.relative_quantize(proposed; bits)
            for k in 1:l
                if iszero(rounded[k, k])
                    rounded[k, k] = only(MS.relative_quantize([proposed[k, k]]; bits))
                end
            end
            rounded
        end
        @assert all(!iszero(T[k, k]) for k in 1:l)
        exactT = Nemo.matrix(Nemo.QQ, T)
        exactM = Nemo.matrix(Nemo.QQ, M)
        C = Nemo.solve(exactT, exactM; side = :right)
        @assert exactT*C == exactM
        result[j] = reshape(nemo_array(C), l, d, di, r)
        previous = result[j-1]
        lp, dp, dip, rp = size(previous)
        @assert rp == l
        result[j-1] = reshape(
            nemo_array(Nemo.matrix(Nemo.QQ, reshape(previous, lp*dp*dip, rp))*exactT),
            lp,
            dp,
            dip,
            l,
        )
        transforms[j] = T
    end
    @assert verify_exact_compression(cores, result, transforms)
    return (; cores = result, transforms)
end

"Certified initial error from three ordinary double-layer contractions."
function initial_error(exact, candidate; precision = 192, weights = WEIGHTS)
    qq = MS.hs_inner_ball(exact, exact; weights, precision)
    cc = MS.hs_inner_ball(candidate, candidate; weights, precision)
    qc = MS.hs_inner_ball(exact, candidate; weights, precision)
    error2 = qq+cc-2qc
    upper2 = max(F(0), PG.relative_upper(error2))
    error = MN.relative_upper_sqrt(upper2)
    @assert error2 < parent(error2)(error^2)
    return (; error, qq, cc, qc, error2)
end

"Exact sparse physical expansion for independent small-system controls."
function sparse_operator(cores)
    partial = [sparse(ones(F, 1, 1))]
    for core in cores
        l, d, _, r = size(core)
        next = [spzeros(F, d*size(first(partial), 1), d*size(first(partial), 2)) for _ in 1:r]
        for a in 1:r, b in 1:l
            C = sparse(core[b, :, :, a])
            (iszero(nnz(C)) || iszero(nnz(partial[b]))) && continue
            next[a] += kron(partial[b], C)
        end
        partial = next
    end
    return only(partial)
end

end
