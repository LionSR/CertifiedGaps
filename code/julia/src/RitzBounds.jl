"""
    CertifiedGaps.RitzBounds

Rigorous upper bounds on the lowest eigenvalues of large sparse symmetric blocks, and from them
rigorous lower bounds on partition functions of rings too large for complete spectra.

By the Courant–Fischer theorem (Cauchy interlacing), the sorted Rayleigh–Ritz values θ₁ ≤ θ₂ ≤ …
of a block M on any k-dimensional subspace S satisfy θᵢ ≥ λᵢ(M), i = 1, …, k. Hence, for any
choice of blocks B with exact multiplicities m_B and of trial subspaces,

    Z(β) = Σ_B m_B Σ_i e^{−βλᵢ(B)} ≥ Σ_B m_B Σ_{i ≤ k_B} e^{−β uᵢ(B)}

for any upper bounds uᵢ(B) ≥ θᵢ(B). Omitting blocks or states, or a poor subspace, only lowers
the bound: the eigensolver (KrylovKit Lanczos) is never trusted. `ritz_upper_bounds` turns an
arbitrary Float64 basis V into rigorous bounds uᵢ with explicit a priori IEEE binary64 rounding
bounds, in the model of `VerifiedSpectra.verified_eigen` (round-to-nearest, gradual underflow,
|fl(Σ aⱼbⱼ) − Σ aⱼbⱼ| ≤ γ_K Σ|aⱼbⱼ| for any summation order and FMA). `partial_ring_spectrum`
applies it to the sparse charge and momentum blocks of `VerifiedSpectra.ring_blocks`, and
`log_partition_lower` sums the result in interval arithmetic. `Seeds.coupled_seed_defects` turns
such lower bounds on log Z_{2n₀} into seed defects.
"""
module RitzBounds

using LinearAlgebra, SparseArrays, IntervalArithmetic, KrylovKit, Random
using ..ExactArithmetic: RationalBig, parse_rational
using ..Results: write_result, read_result
using ..RingModels: RingModel
using ..VerifiedSpectra: verified_eigen, ring_blocks, _gamma, _fnorm, TINY
using ..Recurrence: to_interval

export PartialBlock, ritz_upper_bounds, lowest_upper_bounds, partial_ring_spectrum,
    log_partition_lower, write_partial_spectrum, read_partial_spectrum,
    partial_spectrum_summary

"""
    PartialBlock(upper, multiplicity, label)

`multiplicity` copies of a symmetry block whose i-th smallest eigenvalue is at most `upper[i]`
(sorted, `Float64` rounded up from exact rational bounds) for i ≤ length(upper). The block may
have more eigenvalues; they are not counted.
"""
struct PartialBlock
    upper::Vector{Float64}
    multiplicity::Int
    label::Dict{String,Any}
end

# Smallest Float64 not below the exact rational q.
function _float_up(q::RationalBig)
    f = Float64(q)
    while RationalBig(f) < q
        f = nextfloat(f)
    end
    return f
end

# Row data of a symmetric Float64 matrix: the largest number p of stored entries in a row, the
# computed row sums Σⱼ|Mᵢⱼ|, and the computed Gershgorin lower ends Mᵢᵢ − Σ_{j≠i}|Mᵢⱼ|.
function _rows(M::SparseMatrixCSC{Float64})
    N = size(M, 1)
    p, sums, gersh = 1, zeros(N), zeros(N)
    rows, vals = rowvals(M), nonzeros(M)
    for j in 1:N                                  # column j is row j (M is symmetric)
        range = nzrange(M, j)
        p = max(p, length(range))
        diag, off = 0.0, 0.0
        for k in range
            rows[k] == j ? (diag = vals[k]) : (off += abs(vals[k]))
        end
        sums[j] = abs(diag) + off
        gersh[j] = diag - off
    end
    return p, sums, gersh
end
function _rows(M::Matrix{Float64})
    N = size(M, 1)
    sums = [sum(abs, view(M, :, j)) for j in 1:N]
    gersh = [M[j, j] - (sums[j] - abs(M[j, j])) for j in 1:N]
    return N, sums, gersh
end

"""
    ritz_upper_bounds(M, V; matrix_error = 0) -> NamedTuple

Rigorous upper bounds `upper[i] ≥ λᵢ(A)`, i = 1, …, k, on the k smallest eigenvalues of every
real symmetric A with ‖A − M‖₂ ≤ `matrix_error`, from the Rayleigh–Ritz values of the span of
the columns of an arbitrary `Float64` matrix V (N × k). `M` is an exactly symmetric `Float64`
matrix, sparse or dense. Returns exact rationals `upper` (sorted), their `Float64` round-ups
`upper_float`, the computed Ritz values `ritz` (approximations, not bounds), and the exact
rational ingredients `shift` (c), `orthogonality` (η), `product_error` and `matrix_error`.

Derivation. Write A₀ = VᵀMV and G = VᵀV (exact real numbers, V being the stored binary matrix).

1. Rounding. With p the largest number of stored entries in a row of M and ρ ≥ ‖|M|‖_∞ ≥ ‖|M|‖₂
   (M symmetric), Ŵ = fl(MV) obeys ‖Ŵ − MV‖_F ≤ γ_{p+2} ρ ‖V‖_F, and Â = fl(VᵀŴ) obeys
   ‖Â − VᵀŴ‖_F ≤ γ_{N+2} ‖V‖_F ‖Ŵ‖_F, so ‖Â − A₀‖_F ≤ ε_A := γ_{N+2}‖V‖_F‖Ŵ‖_F +
   ‖V‖_F γ_{p+2} ρ ‖V‖_F (plus underflow allowances). The upper triangle of Â, mirrored, is an
   exactly symmetric S with ‖S − A₀‖_F ≤ √2 ε_A, and `verified_eigen(S; matrix_error = √2 ε_A)`
   gives μᵢ⁺ ≥ μᵢ := λᵢ(A₀). Likewise ‖G − I‖₂ ≤ η := ‖fl(VᵀV) − I‖_F + γ_{N+2}‖V‖_F²; we need
   η < 1/2 (V is then of full rank).
2. Shift. Gershgorin, evaluated with rounding bound γ_{p+2}ρ, gives c ≤ λ_min(M), so
   B := A₀ − cG = Vᵀ(M − c)V ⪰ 0.
3. Pencil. The Ritz values θᵢ of span V are the eigenvalues of the pencil (A₀, G), and
   θᵢ ≥ λᵢ(M) by Courant–Fischer restricted to span V. Then θᵢ − c = λᵢ(G^{−1/2} B G^{−1/2}),
   and by Ostrowski's theorem λᵢ(G^{−1/2} B G^{−1/2}) = tᵢ λᵢ(B) with
   tᵢ ∈ [1/λ_max(G), 1/λ_min(G)] ⊂ [1/(1 + η), 1/(1 − η)]. As λᵢ(B) ≥ 0,
   θᵢ − c ≤ λᵢ(B)/(1 − η), and by Weyl, B = (A₀ − cI) − c(G − I) gives
   λᵢ(B) ≤ μᵢ − c + |c|η. Hence

       λᵢ(A) ≤ λᵢ(M) + ‖A − M‖₂ ≤ θᵢ + matrix_error
             ≤ c + (μᵢ⁺ − c + |c| η)/(1 − η) + matrix_error =: upper[i].

No step trusts the eigensolver that produced V; a poor V only gives weaker bounds.
"""
function ritz_upper_bounds(M::Union{SparseMatrixCSC{Float64},Matrix{Float64}},
    V::Matrix{Float64}; matrix_error = 0)
    N, k = size(V)
    size(M) == (N, N) || throw(DimensionMismatch("M is $(size(M)), V has $N rows"))
    1 <= k <= N || throw(ArgumentError("need 1 ≤ k ≤ N basis vectors, got k = $k"))
    issymmetric(M) || throw(ArgumentError("the block center is not exactly symmetric"))
    all(isfinite, V) || throw(ArgumentError("V has non-finite entries"))
    p, sums, gersh = _rows(M)
    γp, γN = _gamma(p + 2), _gamma(N + 2)
    ρ = RationalBig(maximum(sums)) / (1 - γp) + p * TINY
    c = RationalBig(minimum(gersh)) - γp * ρ - p * TINY
    W = M * V
    Â = transpose(V) * W
    G = transpose(V) * V
    all(isfinite, W) && all(isfinite, Â) && all(isfinite, G) || error("non-finite products")
    # fl(Gᵢᵢ − 1) is exact (Sterbenz) for Gᵢᵢ ∈ [1/2, 2], so D is the exact G − I of fl(VᵀV).
    all(g -> 1 / 2 <= g <= 2, diag(G)) || error("basis vectors far from unit length")
    D = G - I
    nV, nW = _fnorm(V), _fnorm(W)
    εW = γp * ρ * nV + N * k * p * TINY
    εA = γN * nV * nW + nV * εW + k^2 * (N + 2) * TINY
    η = _fnorm(D) + γN * nV^2 + k^2 * (N + 2) * TINY
    η < 1 // 2 || error("basis far from orthonormal (‖VᵀV − I‖ ≤ $(Float64(η)))")
    S = Matrix(Symmetric(Â, :U))
    errS = RationalBig(sup(sqrt(interval(BigFloat, 2)))) * εA
    e = verified_eigen(S; matrix_error = errS)
    ε = RationalBig(matrix_error)
    upper = [c + (RationalBig(μ) + e.radius - c + abs(c) * η) / (1 - η) + ε for μ in e.values]
    issorted(upper) || error("upper bounds not sorted")
    return (; upper, upper_float = _float_up.(upper), ritz = e.values, shift = c,
        orthogonality = η, product_error = errS, matrix_error = ε)
end

"""
    lowest_upper_bounds(M, howmany; matrix_error = 0, dense_max = 1000, tol = 1e-10,
        krylovdim = max(2howmany + 20, 64), maxiter = 200, seed = 1) -> NamedTuple

Rigorous upper bounds on the `howmany` smallest eigenvalues of a block with center `M` (sparse
or dense, exactly symmetric; the true block is within `matrix_error` in spectral norm).
A block of dimension at most `dense_max` is enclosed completely by `verified_eigen` (every
eigenvalue, `howmany` is ignored; complex Hermitian blocks are allowed here). A larger real block
gets approximate eigenvectors from KrylovKit's Lanczos `eigsolve` (with a deterministic start
vector from `seed`), turned into bounds by `ritz_upper_bounds`; non-convergence costs accuracy,
never rigor. Returns `upper` (sorted `Float64` upper bounds), `method` ("dense" or "ritz"), and
for "ritz" the Ritz record (`orthogonality`, `product_error`, `shift`) and the Lanczos residuals
`residual_float` (diagnostics only).
"""
function lowest_upper_bounds(M::AbstractMatrix, howmany::Integer; matrix_error = 0,
    dense_max::Integer = 1000, tol = 1e-10, krylovdim::Integer = max(2howmany + 20, 64),
    maxiter::Integer = 200, seed::Integer = 1)
    N = size(M, 1)
    if N <= dense_max
        e = verified_eigen(Matrix(M); matrix_error)
        upper = _float_up.(RationalBig.(e.values) .+ e.radius)
        return (; upper, method = "dense", radius = e.radius)
    end
    eltype(M) === Float64 ||
        throw(ArgumentError("Ritz bounds for large complex blocks are not implemented"))
    k = min(howmany, N)
    x0 = randn(MersenneTwister(seed), N)
    vals, vecs, info = eigsolve(M, x0, k, :SR; ishermitian = true, krylovdim = min(krylovdim, N),
        tol, maxiter, verbosity = 0)
    V = reduce(hcat, vecs[1:min(k, length(vecs))])
    r = ritz_upper_bounds(M isa SparseMatrixCSC ? M : Matrix{Float64}(M), V; matrix_error)
    return (; upper = r.upper_float, method = "ritz", r.orthogonality, r.product_error, r.shift,
        residual_float = maximum(info.normres[1:min(k, length(info.normres))]),
        converged = info.converged)
end

"""
    partial_ring_spectrum(model, L; howmany, charges = nothing, dense_max = 1000,
        progress = false, kwargs...) -> Vector{PartialBlock}

Rigorous upper bounds on the low eigenvalues of every block of the ring of L sites in the
sectors `charges` (all if `nothing`), from the sparse exact blocks of
`VerifiedSpectra.ring_blocks(model, L; sparse = true, charges)` (interval-enclosed entries,
exact `matrix_error`) and `lowest_upper_bounds` (`kwargs` are passed on). `howmany` is an
integer or a function of the block label returning the number of levels to bound (0 skips the
block). Each block keeps its exact multiplicity, so `log_partition_lower` gives a rigorous lower
bound on log Z_L(β).
"""
function partial_ring_spectrum(model::RingModel, L::Integer; howmany, charges = nothing,
    dense_max::Integer = 1000, progress = false, kwargs...)
    out = PartialBlock[]
    count(label) = howmany isa Integer ? howmany : howmany(label)
    ring_blocks(model, L; sparse = true, charges) do M, ε, mult, label
        k = count(label)
        k > 0 || return
        seconds = @elapsed r = lowest_upper_bounds(M, k; matrix_error = ε, dense_max, kwargs...)
        label["seconds"] = seconds
        label["method"] = r.method
        label["matrix_error"] = string(ε)
        if r.method == "ritz"
            label["orthogonality"] = string(r.orthogonality)
            label["product_error"] = string(r.product_error)
            label["residual_float"] = r.residual_float
            label["converged"] = r.converged
        end
        push!(out, PartialBlock(r.upper, mult, label))
        progress && println("L=$L charge=$(label["charge"]) m=$(label["momentum"]) ",
            "dim=$(size(M, 1)) levels=$(length(r.upper)) $(r.method) ",
            "low=$(r.upper[1]) $(round(seconds; digits = 1))s")
        flush(stdout)
    end
    return out
end

"""
    log_partition_lower(blocks, β) -> Interval{Float64}

Enclosure of log Σ_B m_B Σᵢ e^{−β upperᵢ(B)} for exact or interval β ≥ 0. Its lower end is a
rigorous lower bound on log Z(β) of the ring; its upper end bounds only the truncated sum.
"""
function log_partition_lower(blocks::AbstractVector{PartialBlock}, β)
    b = β isa Interval ? β : to_interval(RationalBig(β))
    inf(b) >= 0 || throw(ArgumentError("β must be nonnegative"))
    shift = minimum(minimum(s.upper) for s in blocks if !isempty(s.upper))
    s0 = interval(shift)
    total = interval(0.0)
    for s in blocks, u in s.upper
        total += s.multiplicity * exp(-b * (interval(u) - s0))
    end
    return log(total) - b * s0
end

const PARTIAL_FORMAT = "CertifiedGaps partial spectrum upper bounds v1"

"""
    write_partial_spectrum(path, blocks) -> path

Store the upper bounds losslessly (shortest round-trip decimals) with multiplicities and labels,
via `Results.write_result` (`.json` or `.json.gz`); the file is read back and compared.
"""
function write_partial_spectrum(path, blocks::AbstractVector{PartialBlock})
    data = [Dict{String,Any}("upper" => string.(s.upper), "multiplicity" => s.multiplicity,
        "label" => s.label) for s in blocks]
    write_result(path, Dict{String,Any}("format" => PARTIAL_FORMAT, "blocks" => data);
        compact = true)
    back = read_partial_spectrum(path)
    length(back) == length(blocks) && all(a.upper == b.upper && a.multiplicity ==
        b.multiplicity for (a, b) in zip(back, blocks)) ||
        error("$path does not read back as the bounds written")
    return path
end

"Inverse of `write_partial_spectrum`."
function read_partial_spectrum(path)
    data = read_result(path; exact = false)
    data["format"] == PARTIAL_FORMAT || error("$path is not a partial spectrum file")
    return [PartialBlock(parse.(Float64, b["upper"]), Int(b["multiplicity"]),
        Dict{String,Any}(String(k) => v for (k, v) in b["label"])) for b in data["blocks"]]
end

"JSON-ready summary of a partial spectrum (counts, dimensions covered, times)."
function partial_spectrum_summary(blocks::AbstractVector{PartialBlock})
    Dict{String,Any}(
        "blocks" => length(blocks),
        "levels" => sum(length(s.upper) * s.multiplicity for s in blocks),
        "dimension_of_blocks" => sum(s.label["dimension"] * s.multiplicity for s in blocks),
        "ritz_blocks" => count(s -> s.label["method"] == "ritz", blocks),
        "lowest_upper_bound" => minimum(minimum(s.upper) for s in blocks),
        "seconds" => sum(s.label["seconds"] for s in blocks),
    )
end

end
