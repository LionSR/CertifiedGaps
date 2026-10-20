"""
    TrialSeed

Untrusted generator of the data behind the kernel-checked trial-state bound
`lean/SeedKernel/Heisenberg/Trial*.lean`: the integer matrix-product trial state `A` of
`models/spin1/src/trials.jl`, its transfer matrix `T = Σₛ Aₛ ⊗ Aₛ` and energy insertion
`O = Σ h(st; s't') (AₛAₜ) ⊗ (Aₛ'Aₜ')` on the 676-dimensional doubled virtual space, and an integer
intertwiner `W` with `T W = W B`, `O W = W C` and `V W = c 1`, where `B` and `C` are block diagonal
with one block `1 ⊗ tᵦ` (resp. `1 ⊗ oᵦ`) per spin-`J` sector and sign of the swap symmetry.

Nothing computed here is trusted: Lean rechecks every identity. The module only has to produce
data that passes those checks.

Conventions match `CertifiedGaps.SpinOneHeisenberg` (`TrialBound.lean`, `SzBridge.lean`): the label
`s ∈ -1:1` is the index `s + 1 : Fin 3`, the virtual labels keep the order of `trial_matrices()`,
and the doubled index `(a, b) : Fin 26 × Fin 26` is the integer `26a + b` (0-based).
"""
module TrialSeed

using LinearAlgebra, SparseArrays
import Nemo
using CertifiedGaps.Results: read_exact_json

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SPIN1_SRC = joinpath(REPO, "models", "spin1", "src")
"Upstream evidence directory read by `trial_matrices()` (fetched, not committed)."
const SOURCE = normpath(joinpath(REPO, "references", "periodic-gap-source", "verification",
    "computations"))
integer(x::AbstractString) = parse(BigInt, x)
integer(x::Integer) = BigInt(x)
integer(x) = error("Refusing inexact integer input: $(typeof(x))")
# `trial_matrices` (the 23-integer recipe, asserted equal to upstream `matrix-data.json`) and
# `virtual_ladders` (the virtual spin ladders and their commutators with `A`).
include(joinpath(SPIN1_SRC, "trials.jl"))
include(joinpath(SPIN1_SRC, "trials_su2.jl"))

const D = 26

"Sparse integer Kronecker product with the `(a, b) ↦ 26a + b` doubled index."
kr(X, Y) = kron(sparse(X), sparse(Y))

"""
    operators()

The recipe matrices `A[s]` (`s ∈ -1:1`), the virtual ladders, and the doubled-space operators
`T`, `O`, raising `R = 2 plus ⊗ 1 - 1 ⊗ minus`, lowering `L = minus ⊗ 1 - 2 (1 ⊗ plus)` and the
swap `SW`, all as sparse `Int` matrices (entries are far below `2^62`).
"""
function operators()
    labels, Abig = trial_matrices()
    plus, minus = virtual_ladders(labels, Abig)
    A = Dict(s => Int.(Abig[s]) for s in -1:1)
    one26 = Matrix{Int}(I, D, D)
    T = sum(kr(A[s], A[s]) for s in -1:1)
    X = Dict((s, t) => A[s] * A[t] for s in -1:1, t in -1:1)
    # h = SᶻSᶻ + (S⁺S⁻ + S⁻S⁺)/2 in the Sᶻ basis, an integer matrix.
    O = sum(s * t * kr(X[(s, t)], X[(s, t)]) for s in -1:1, t in -1:1)
    for s in 0:1, t in -1:0
        O += kr(X[(s, t)], X[(s - 1, t + 1)]) + kr(X[(s - 1, t + 1)], X[(s, t)])
    end
    R = kr(2Int.(plus), one26) - kr(one26, Int.(minus))
    L = kr(Int.(minus), one26) - 2kr(one26, Int.(plus))
    SW = spzeros(Int, D^2, D^2)
    for a in 0:D-1, b in 0:D-1
        SW[b*D+a+1, a*D+b+1] = 1
    end
    @assert T * L == L * T && T * R == R * T && O * L == L * O && O * R == R * O
    @assert SW * T * SW == T && SW * O * SW == O && SW * R * SW == -L
    return (; labels, A, T, O, R, L, SW)
end

qq(M) = Nemo.matrix(Nemo.QQ, Matrix{BigInt}(M))
tobig(x) = BigInt(Nemo.numerator(x)) // BigInt(Nemo.denominator(x))

"Smallest positive integer multiple of a rational matrix (columns scaled together)."
function integral(M::Nemo.QQMatrix)
    d = foldl(lcm, (BigInt(Nemo.denominator(M[i, j])) for i in 1:Nemo.nrows(M),
        j in 1:Nemo.ncols(M)); init=big(1))
    Z = [BigInt(Nemo.numerator(M[i, j] * d)) for i in 1:Nemo.nrows(M), j in 1:Nemo.ncols(M)]
    g = foldl(gcd, Z; init=big(0))
    return iszero(g) ? Z : Z .÷ g
end

"Solve `N x = Y` exactly for a full-column-rank `N`; asserts the solution is integral."
function solve_integral(N::Matrix{BigInt}, Y::Matrix{BigInt})
    x = Nemo.solve(qq(N), qq(Y); side=:right)
    @assert qq(N) * x == qq(Y)
    out = [tobig(x[i, j]) for i in 1:Nemo.nrows(x), j in 1:Nemo.ncols(x)]
    @assert all(isone(denominator(v)) for v in out) "non-integral block"
    return numerator.(out)
end

"Integer basis of the saturated lattice `ker M ∩ ℤⁿ` for an integer matrix `M`."
function lattice_kernel(M::Matrix{BigInt})
    Z = Nemo.matrix(Nemo.ZZ, M)
    K = Nemo.kernel(Z; side=:right)
    # `kernel` over ℤ returns a basis of the kernel as a ℤ-module (via HNF), hence saturated.
    return [BigInt(K[i, j]) for i in 1:Nemo.nrows(K), j in 1:Nemo.ncols(K)]
end

"""
    blocks(ops)

Highest-weight bases split by the swap symmetry. For each spin `J` the highest-weight vectors
`N` (weight `2J`, killed by `R`) carry the integer block `t` (`T N = N t`) and the involution
`M` (`R^{2J} SW N = N M`, `M² = λ² 1`) commuting with `t`; the integer kernel lattices of
`M ∓ λ` give the blocks `(J, ±)`, each with integer `tᵦ`, `oᵦ`.
"""
function blocks(ops)
    (; labels, T, O, R, SW) = ops
    q(i) = labels[div(i - 1, D)+1][3] - labels[mod(i - 1, D)+1][3]
    out = Any[]
    for J in 0:5
        sec = [i for i in 1:D^2 if q(i) == 2J]
        N = zeros(BigInt, D^2, 0)
        if J == 5
            N = zeros(BigInt, D^2, 1)
            N[sec[1], 1] = 1
        else
            _, K = Nemo.nullspace(qq(R[:, sec]))
            Kz = integral(K)
            N = zeros(BigInt, D^2, size(Kz, 2))
            N[sec, :] = Kz
        end
        X = BigInt.(SW) * N
        for _ in 1:2J
            X = BigInt.(R) * X
        end
        M = solve_integral(N, X)
        λ2 = (M*M)[1, 1]
        @assert M * M == λ2 * I
        λ = isqrt(λ2)
        @assert λ^2 == λ2
        for sign in (1, -1)
            K = lattice_kernel(M - sign * λ * I)
            size(K, 2) == 0 && continue
            H = N * K
            t = solve_integral(H, BigInt.(T) * H)
            o = solve_integral(H, BigInt.(O) * H)
            push!(out, (; J, sign, H, t, o))
        end
    end
    @assert sum((2b.J + 1) * size(b.t, 1) for b in out) == D^2
    return out
end

"""
    intertwiner(ops, bs)

The columns `L^k Hᵦ[:, c] / gᵦₖ` in the order `(b, k, c)` (block, lowering power, copy), each
`(b, k)` group divided by its gcd; the local column index inside block `b` is `k mᵦ + c`.
Returns `W` (sparse), the block offsets, and `V = c W⁻¹` with `c` the least common denominator.
"""
function intertwiner(ops, bs)
    Lb = BigInt.(ops.L)
    cols = Vector{Vector{BigInt}}()
    offsets = Int[]
    for b in bs
        push!(offsets, length(cols))
        X = b.H
        for k in 0:2b.J
            g = foldl(gcd, X; init=big(0))
            Xs = X .÷ g
            append!(cols, [Xs[:, c] for c in 1:size(X, 2)])
            X = Lb * X
        end
    end
    W = sparse(reduce(hcat, cols))
    # W is block diagonal in the charge q(row) = weight(column): invert sector by sector.
    labels = ops.labels
    q(i) = labels[div(i - 1, D)+1][3] - labels[mod(i - 1, D)+1][3]
    weight = Int[]
    for b in bs, k in 0:2b.J, _ in 1:size(b.t, 1)
        push!(weight, 2b.J - 2k)
    end
    Vq = Dict{Tuple{Int,Int},Rational{BigInt}}()
    for δ in -10:2:10
        rows = [i for i in 1:D^2 if q(i) == δ]
        cs = [j for j in 1:D^2 if weight[j] == δ]
        @assert length(rows) == length(cs)
        Winv = inv(qq(Matrix(W[rows, cs])))
        for (a, j) in enumerate(cs), (bidx, i) in enumerate(rows)
            v = Winv[a, bidx]
            iszero(v) || (Vq[(j, i)] = tobig(v))
        end
    end
    c = foldl(lcm, (denominator(v) for v in values(Vq)); init=big(1))
    V = spzeros(BigInt, D^2, D^2)
    for ((j, i), v) in Vq
        V[j, i] = numerator(v * c)
    end
    return (; W, V, c, offsets)
end

"Block-diagonal `B` (or `C`) in the column order of `intertwiner`."
function block_matrix(bs, which)
    n = sum((2b.J + 1) * size(b.t, 1) for b in bs)
    Bm = spzeros(BigInt, n, n)
    off = 0
    for b in bs
        m = size(b.t, 1)
        blk = which == :t ? b.t : b.o
        for k in 0:2b.J
            Bm[off+k*m+1:off+(k+1)*m, off+k*m+1:off+(k+1)*m] = blk
        end
        off += (2b.J + 1) * m
    end
    return Bm
end

"Exact `tr(t^N)` and `tr(t^(N-2) o)` of every block."
function block_traces(bs, N)
    [(V=tr(b.t^N), U=tr(b.t^(N - 2) * b.o)) for b in bs]
end

end
