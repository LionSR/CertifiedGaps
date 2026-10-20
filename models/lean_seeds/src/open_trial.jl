"""
    OpenTrialSeed

Untrusted generator of the data behind the kernel-checked open-chain trial bound
`lean/SeedKernel/Heisenberg/OpenTrial*.lean`: the opposite-field trial vectors of OpenAI's
*A boundary-field gap for the spin-one Heisenberg chain* (appendix "Exact trial vectors for
opposite fields", Proposition `apptr:bound`).

The trial vector on `m` sites is the open matrix-product state `ψ(σ) = bᵀ A_{σ₁} ⋯ A_{σₘ} c` with
the bond-dimension-26 tensors `A` of the periodic trial (`trial.jl`) and the endpoint vectors
`b`, `c` of the paper's endpoint table. Its norm and its energy in

    H_m^opp = ∑ⱼ Sⱼ · Sⱼ₊₁ - (3/5) S₁ᶻ + (3/5) Sₘᶻ

are contractions `L₀ ⋯ u` of the transfer matrix `T`, the bond insertion `O` and the field
insertion `G = ∑ₛ s Aₛ ⊗ Aₛ` between `L₀ = b ⊗ b` and `u = c ⊗ c`. Lean reduces them with the
periodic trial's intertwiner (`T W = W B`, `O W = W C`, `V W = c 1`) to the blocks `tᵦ`, `oᵦ`:
`L₀ W` is supported on the weight-zero copy `k = Jᵦ` of each block, so each block contributes
one row iteration of size `mᵦ ≤ 32`. This module runs exactly the iteration that the kernel runs
(`OpenTrialCompute.lean`), including its lane-width checks, and records the per-block outputs.

Nothing computed here is trusted: Lean recomputes every output. The module only has to produce
data that pass those checks; it also cross-checks the totals against the paper's
`trial-exact-results.json`.
"""
module OpenTrialSeed

include(joinpath(@__DIR__, "trial.jl"))
using .TrialSeed
using SparseArrays, LinearAlgebra
import JSON

const REPO = TrialSeed.REPO
"Upstream evidence of the boundary-field paper (fetched by `tools/fetch_upstream.jl`)."
const BOUNDARY_SOURCE = joinpath(REPO, "references", "boundary-field-gap-source", "verification")
const D = TrialSeed.D

"The lengths `m` of Proposition `apptr:bound`."
const LENGTHS = (34, 78, 142, 240)
"Endpoint parameters `βᵢ`, `κᵢ` of the paper's endpoint table, one per multiplet."
const BETA = [10000, 5620, 2456, 4868, 82, -809, 12475, -9645]
const KAPPA = [12, 12, 12, 12, 6, 6, 6, 1]
"Scaled coefficients: `500000 h_bond`, `500000 · 3/5`, `500000 a`, and `500000 α`."
const BOND, FIELD, SHIFT, ALPHA = 500000, 300000, 700741, 202100

"""
    endpoints(labels)

The endpoint columns `b`, `c` (`b = βᵢ`, `c = κᵢ βᵢ` on the label `d = -1` of multiplet `i`).
"""
function endpoints(labels)
    multiplets = unique([(l, a) for (l, a, _) in labels])
    b = zeros(BigInt, D)
    c = zeros(BigInt, D)
    for (i, (l, a, d)) in enumerate(labels)
        d == -1 || continue
        k = findfirst(==((l, a)), multiplets)
        b[i] = BETA[k]
        c[i] = KAPPA[k] * BETA[k]
    end
    return b, c
end

"The paper's verification data (tensors, endpoints, exact totals)."
function upstream()
    mats = JSON.parsefile(joinpath(BOUNDARY_SOURCE, "trial-exact-matrices.json"))
    res = JSON.parsefile(joinpath(BOUNDARY_SOURCE, "trial-exact-results.json"))
    return mats, res
end

big_int(x::AbstractString) = parse(BigInt, x)
big_int(x::Integer) = BigInt(x)

"""
    direct_totals(ops, b, c)

`(V_m, H_m^*)` by the paper's three-row recurrence on the full doubled space (676 coordinates).
"""
function direct_totals(ops, b, c)
    T = BigInt.(ops.T)
    O = BigInt.(ops.O)
    G = sum(s * TrialSeed.kr(ops.A[s], ops.A[s]) for s in -1:1)
    G = BigInt.(G)
    L0 = permutedims(kron(b, b))
    u = kron(c, c)
    x = L0
    y = L0 * T
    w = -FIELD * (L0 * G)
    out = Dict{Int,Tuple{BigInt,BigInt}}()
    Gu = G * u
    for m in 2:maximum(LENGTHS)
        w, x, y = w * T + BOND * (x * O) + SHIFT * (y * T), y, y * T
        if m in LENGTHS
            V = (y*u)[1]
            H = (w*u)[1] + FIELD * (x*Gu)[1] + SHIFT * V
            out[m] = (V, H)
        end
    end
    return out
end

"Largest absolute entry (zero for an empty list)."
maxabs(xs) = maximum(abs, xs; init=big(0))

"""
    vecmat_ok(W, n, β, off)

The lane condition of the packed vector-matrix product of `OpenTrialCompute.vecmat`:
`n² (2β) (2 off) < 2^W`.
"""
vecmat_ok(W, n, β, off) = n^2 * (2β) * (2off) < big(2)^W

"""
    block_run(t, o, ℓ, ℓ′, ν, ν′; Wt, Wo, offt, offo)

The kernel's row iteration on one block (`OpenTrialCompute.blockRun`): `x_r = ℓ t^r`,
`z₁ = -300000 ℓ′`, `z_{r+1} = z_r t + 500000 x_{r-1} o`; at each length `m` the outputs
`(x_m ν, z_m ν, x_{m-1} ν′)`. Also returns the largest `β` met by each product, from which the
lane widths are chosen, and whether all lane conditions hold for the given widths.
"""
function block_run(t, o, ℓ, ℓ′, ν, ν′; Wt=0, Wo=0, offt=0, offo=0)
    n = length(ℓ)
    βt, βo = big(0), big(0)
    # The lane condition of every product the kernel checks: the input of each product by `t`
    # (`okT`) and by `o` (`okO`), including those of the last state, whose step is never used.
    okT(x) = (β = maxabs(x); βt = max(βt, β); vecmat_ok(Wt, n, β, offt))
    okO(x) = (β = maxabs(x); βo = max(βo, β); vecmat_ok(Wo, n, β, offo))
    vt(x) = permutedims(permutedims(x) * t)[:, 1]
    vo(x) = permutedims(permutedims(x) * o)[:, 1]
    ok = okT(ℓ)
    xp, x, z = ℓ, vt(ℓ), -FIELD * ℓ′
    outs = Tuple{BigInt,BigInt,BigInt}[]
    for r in 1:maximum(LENGTHS)
        r in LENGTHS && push!(outs, (dot(x, ν), dot(z, ν), dot(xp, ν′)))
        ok = okT(x) & okT(z) & okO(xp) & ok
        r == maximum(LENGTHS) && break
        xp, x, z = x, vt(x), vt(z) + BOND * vo(xp)
    end
    return (; outs, βt, βo, ok)
end

"Bit length of a nonnegative integer."
bitlen(x) = x == 0 ? 0 : ndigits(x, base=2)

"""
    reduce_blocks()

Everything the generator writes: per block the lane widths and outputs, and the totals. Asserts
the reduction (zero copies, intertwiner sums) and the agreement with the direct recurrence and
the paper's `trial-exact-results.json`.
"""
function reduce_blocks()
    ops = TrialSeed.operators()
    bs = TrialSeed.blocks(ops)
    tw = TrialSeed.intertwiner(ops, bs)
    W, V, cs = BigInt.(tw.W), tw.V, tw.c
    labels = ops.labels
    b, c = endpoints(labels)

    # The tensors and endpoints are the paper's.
    mats, res = upstream()
    for s in -1:1, i in 1:D, j in 1:D
        @assert ops.A[s][i, j] == big_int(mats["A"][string(s)][i][j])
    end
    @assert b == big_int.(mats["b"]) && c == big_int.(mats["c"])

    G = BigInt.(sum(s * TrialSeed.kr(ops.A[s], ops.A[s]) for s in -1:1))
    L0 = kron(b, b)
    u = kron(c, c)
    xW = permutedims(permutedims(L0) * W)[:, 1]
    xGW = permutedims(permutedims(L0) * G * W)[:, 1]
    Vu = V * u
    VGu = V * (G * u)

    direct = direct_totals(ops, b, c)
    for m in LENGTHS
        @assert direct[m][1] == big_int(res[string(m)]["V_m"])
        @assert direct[m][2] == big_int(res[string(m)]["H_m_star"])
    end

    blocks = Any[]
    totV = Dict(m => big(0) for m in LENGTHS)
    totH = Dict(m => big(0) for m in LENGTHS)
    for (bi, blk) in enumerate(bs)
        n = size(blk.t, 1)
        off = tw.offsets[bi]
        J = blk.J
        seg(k) = off+k*n+1:off+(k+1)*n
        # Only the weight-zero copy `k = J` of `L₀ W` and `L₀ G W` is nonzero.
        for k in 0:2J
            k == J && continue
            @assert all(iszero, xW[seg(k)]) && all(iszero, xGW[seg(k)])
        end
        ℓ, ℓ′ = xW[seg(J)], xGW[seg(J)]
        ν, ν′ = Vu[seg(J)], VGu[seg(J)]
        t, o = BigInt.(blk.t), BigInt.(blk.o)
        offt, offo = maxabs(t), maxabs(o)
        probe = block_run(t, o, ℓ, ℓ′, ν, ν′)
        Wt = bitlen(n^2 * 4 * probe.βt * offt) + 1
        Wo = bitlen(n^2 * 4 * probe.βo * offo) + 1
        run = block_run(t, o, ℓ, ℓ′, ν, ν′; Wt, Wo, offt, offo)
        @assert run.ok
        for (m, (p1, q, p4)) in zip(LENGTHS, run.outs)
            totV[m] += p1
            totH[m] += q + FIELD * p4 + SHIFT * m * p1
        end
        push!(blocks, (; J, n, Wt, Wo, outs=run.outs))
    end
    for m in LENGTHS
        @assert totV[m] == cs * direct[m][1]
        @assert totH[m] == cs * direct[m][2]
        @assert totV[m] > 0 && totH[m] < ALPHA * totV[m]
    end
    return (; blocks, totals=direct, scale=cs)
end

end
