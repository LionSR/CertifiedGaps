"""
    CertifiedGaps.VerifiedSpectra

Rigorous enclosures of the complete spectrum of a `RingModel` on a ring of L sites, and of
log Z_L(β) from them.

The ring is split exactly into charge sectors and lattice-momentum blocks. Each block matrix is
known as an interval matrix (exact rational couplings times enclosed phases e^{2πiml/L} and
norm ratios); its floating-point center M is diagonalised by LAPACK, and the result is turned
into a guarantee:

* with X the computed eigenvectors and Λ the computed eigenvalues (both exact binary numbers),
  R = MX − XΛ and G = X†X − I are evaluated in Float64 with an a priori rounding-error bound
  (|fl(Σ aᵢbᵢ) − Σ aᵢbᵢ| ≤ γ_K Σ|aᵢbᵢ|, valid for any summation order and with FMA);
* Weyl: B = X†MX = Λ + GΛ + X†R, so |λᵢ(B) − Λᵢ| ≤ ‖G‖‖Λ‖ + ‖X‖‖R‖ (sorted eigenvalues);
* Ostrowski: λᵢ(B) = θᵢ λᵢ(M) with θᵢ ∈ [1 − ‖G‖, 1 + ‖G‖];
* Weyl again for the distance between M and the true block.

Hence every true eigenvalue lies within one radius r of the sorted computed values. No
inequality rests on an unverified floating-point result. Scalar bounds are exact rationals.
"""
module VerifiedSpectra

using LinearAlgebra, IntervalArithmetic
using ..ExactArithmetic: RationalBig
using ..RingModels: RingModel, isreal_model, _reduce
using ..Recurrence: to_interval

export BlockSpectrum, enclose_hermitian, ring_blocks, verified_ring_spectrum,
    log_partition_enclosure, spectrum_summary

const UNIT = big(1) // big(2)^53                 # unit roundoff of binary64
const TINY = big(1) // big(2)^1021              # covers gradual underflow or flush-to-zero

"""
    BlockSpectrum(centers, radius, multiplicity, label)

`multiplicity` copies of a spectrum whose sorted eigenvalues lie within `radius` (an exact
rational) of the sorted `centers`.
"""
struct BlockSpectrum
    centers::Vector{Float64}
    radius::RationalBig
    multiplicity::Int
    label::Dict{String,Any}
end

_gamma(k) = (@assert k * UNIT < 1 // 2; k * UNIT / (1 - k * UNIT))

# Upper bound on the Frobenius norm of a Float64 array, as an exact rational.
function _fnorm(A::AbstractArray{Float64})
    n = length(A)
    s = sum(abs2, A; init = 0.0)
    @assert isfinite(s)
    return _sqrt_up(RationalBig(s) / (1 - _gamma(n + 1)) + n * TINY)
end
_sqrt_up(q::RationalBig) = RationalBig(sup(sqrt(interval(BigFloat, q))))

"""
    enclose_hermitian(M; matrix_error = 0) -> (centers, radius)

Sorted eigenvalue approximations of the Hermitian Float64 matrix `M` (real or complex) and an
exact rational `radius` such that the sorted eigenvalues of every Hermitian A with
‖A − M‖₂ ≤ `matrix_error` lie within `radius` of `centers`, index by index.
"""
function enclose_hermitian(M::Matrix{T}; matrix_error = 0) where {T<:Union{Float64,ComplexF64}}
    N = size(M, 1)
    F = eigen(Hermitian(M))
    Λ, X = F.values, F.vectors
    @assert all(isfinite, Λ) && all(isfinite, X)
    K = N + 8                                     # roundings along any dependency path
    if T === Float64
        R = M * X - X .* transpose(Λ)
        G = transpose(X) * X - I
        nM, nX = _fnorm(M), _fnorm(X)
        λmax = RationalBig(maximum(abs, Λ))
        r = _fnorm(R) + _gamma(K) * (nM * nX + nX * λmax) + N * K * TINY
        η = _fnorm(G) + _gamma(K) * (nX^2 + N) + N * K * TINY
    else
        Mr, Mi, Xr, Xi = real(M), imag(M), real(X), imag(X)
        # Four real GEMMs per product; no complex BLAS kernel is involved.
        Rr = Mr * Xr - Mi * Xi - Xr .* transpose(Λ)
        Ri = Mr * Xi + Mi * Xr - Xi .* transpose(Λ)
        Gr = transpose(Xr) * Xr + transpose(Xi) * Xi - I
        Gi = transpose(Xr) * Xi - transpose(Xi) * Xr
        nMr, nMi, nXr, nXi = _fnorm(Mr), _fnorm(Mi), _fnorm(Xr), _fnorm(Xi)
        λmax = RationalBig(maximum(abs, Λ))
        γ = _gamma(2N + 8)
        r = _fnorm(Rr) + _fnorm(Ri) +
            γ * (nMr * nXr + nMi * nXi + nXr * λmax + nMr * nXi + nMi * nXr + nXi * λmax) +
            2N * (2N + 8) * TINY
        η = _fnorm(Gr) + _fnorm(Gi) + γ * (2 * (nXr^2 + nXi^2) + N) + 2N * (2N + 8) * TINY
    end
    η < 1 // 4 || error("eigenvectors far from orthonormal (‖X†X − I‖ ≤ $(Float64(η)))")
    # |λᵢ(M) − Λᵢ| ≤ (‖G‖ max|Λ| + ‖X‖ ‖R‖ + ‖G‖ |Λᵢ|)/(1 − ‖G‖), with ‖X‖ ≤ 1 + ‖G‖.
    radius = (2η * λmax + (1 + η) * r) / (1 - η) + RationalBig(matrix_error)
    return sort(Λ), radius
end

# State codes: digit of site x (0-based) is (code ÷ d^x) % d. The translation T moves the
# digit at site x to site x+1 (mod L).
_translate(code, d, top) = (code % top) * d + code ÷ top

"""
    ring_blocks(f, model, L)

Exact charge-sector and momentum decomposition of the ring of L sites. For each block it calls
`f(M, matrix_error, multiplicity, label)` with the Float64 center `M` of the interval block
matrix and an exact bound `matrix_error` on ‖A − M‖₂ for the true block A. Blocks are built
one at a time, so only one is in memory. For a real model the momenta ±k have the same
spectrum (A₋ₖ = conj Aₖ), so only 0 ≤ m ≤ L/2 is built, with multiplicity 2 for 0 < m < L/2.
"""
function ring_blocks(f, model::RingModel, L::Integer)
    L >= 2 || throw(ArgumentError("a ring needs L ≥ 2"))
    d = model.d
    dim = big(d)^L
    dim <= 2^26 || throw(ArgumentError("ring dimension $dim is too large for this method"))
    n, top = Int(dim), d^(L - 1)
    codes_by_charge = Dict{Int,Vector{Int}}()
    for code in 0:(n-1)
        q, c = 0, code
        for _ in 1:L
            q += model.charge[c%d+1]
            c ÷= d
        end
        push!(get!(codes_by_charge, _reduce(model, q), Int[]), code)
    end
    real_model = isreal_model(model)
    momenta = real_model ? (0:(L÷2)) : (0:(L-1))
    total = 0
    for q in sort!(collect(keys(codes_by_charge)))
        table = _sector_table(model, L, codes_by_charge[q], n, top)
        for m in momenta
            built = _momentum_block(table, L, m, real_model)
            built === nothing && continue
            M, ε = built
            mult = real_model && 0 < 2m < L ? 2 : 1
            total += size(M, 1) * mult
            f(M, ε, mult, Dict{String,Any}("charge" => q, "momentum" => m,
                "dimension" => size(M, 1)))
        end
    end
    total == n || error("blocks cover $total of $n states")
    return nothing
end

# Orbit representatives of one sector and the exact action of H on them: for representative
# r, the coefficients c with H r = Σ c · T^l s over representatives s and shifts l.
function _sector_table(model::RingModel, L, codes, n, top)
    d = model.d
    rep = Dict{Int,Tuple{Int,Int}}()             # code => (representative index, shift l)
    reps, periods = Int[], Int[]
    for code in codes
        haskey(rep, code) && continue
        orbit = [code]
        while (next = _translate(orbit[end], d, top)) != code
            push!(orbit, next)
        end
        t = argmin(orbit) - 1                    # T^t code = representative
        p = length(orbit)
        push!(reps, orbit[t+1])
        push!(periods, p)
        for (j, s) in enumerate(orbit)
            rep[s] = (length(reps), mod(j - 1 - t, p))
        end
    end
    T = isreal_model(model) ? Rational{Int} : Complex{Rational{Int}}
    convert_entry(x) = T === Rational{Int} ? Rational{Int}(real(x)) :
                       Complex{Rational{Int}}(Rational{Int}(real(x)), Rational{Int}(imag(x)))
    bond_cols = [[(i, convert_entry(model.bond[i, j])) for i in 1:d^2 if !iszero(model.bond[i, j])]
                 for j in 1:d^2]
    site_cols = [[(i, convert_entry(model.onsite[i, j])) for i in 1:d if !iszero(model.onsite[i, j])]
                 for j in 1:d]
    powers = [d^x for x in 0:(L-1)]
    columns = Vector{Dict{Tuple{Int,Int},T}}(undef, length(reps))
    for (k, r) in enumerate(reps)
        column = Dict{Tuple{Int,Int},T}()
        digit(x) = (r ÷ powers[x+1]) % d
        for x in 0:(L-1)
            y = (x + 1) % L
            a, b = digit(x), digit(y)
            for (i, c) in bond_cols[a*d+b+1]
                a2, b2 = divrem(i - 1, d)
                s = r + (a2 - a) * powers[x+1] + (b2 - b) * powers[y+1]
                key = rep[s]
                column[key] = get(column, key, zero(T)) + c
            end
            for (i, c) in site_cols[a+1]
                s = r + (i - 1 - a) * powers[x+1]
                key = rep[s]
                column[key] = get(column, key, zero(T)) + c
            end
        end
        columns[k] = column
    end
    return (; periods, columns)
end

const _CI = Complex{Interval{Float64}}

# Interval block at momentum k = 2πm/L, then its center and an exact bound on ‖A − M‖_F:
#     ⟨s,k|H|r,k⟩ = √(p_r/p_s) Σ_l c(r → T^l s) e^{ikl}.
function _momentum_block(table, L, m, real_model)
    allowed = findall(p -> (m * p) % L == 0, table.periods)
    isempty(allowed) && return nothing
    index = zeros(Int, length(table.periods))
    index[allowed] = 1:length(allowed)
    N = length(allowed)
    realblock = real_model && (2m) % L == 0
    M = realblock ? zeros(Float64, N, N) : zeros(ComplexF64, N, N)
    phases = [complex(cospi(interval(Float64, 2 * ((m * l) % L) // L)),
        sinpi(interval(Float64, 2 * ((m * l) % L) // L))) for l in 0:(L-1)]
    ratios = Dict{Tuple{Int,Int},Interval{Float64}}()
    coefficient = Dict{Any,_CI}()
    enclose(c) = get!(coefficient, c) do
        complex(to_interval(Rational{BigInt}(real(c))), to_interval(Rational{BigInt}(imag(c))))
    end
    err2 = interval(0.0)
    for (j, r) in enumerate(allowed)
        entries = Dict{Int,_CI}()
        for ((s, l), c) in table.columns[r]
            i = index[s]
            (i == 0 || i > j) && continue          # upper triangle, allowed states only
            ratio = get!(ratios, (table.periods[r], table.periods[s])) do
                sqrt(interval(Float64, table.periods[r]) / interval(Float64, table.periods[s]))
            end
            entries[i] = get(entries, i, complex(interval(0.0), interval(0.0))) +
                         enclose(c) * phases[l+1] * ratio
        end
        for (i, v) in entries
            re = mid(real(v))
            im = (realblock || i == j) ? 0.0 : mid(imag(v))
            if realblock
                M[i, j] = re
                M[j, i] = re
            else
                M[i, j] = complex(re, im)
                M[j, i] = complex(re, -im)
            end
            e = interval(mag(real(v) - interval(re)))^2 + interval(mag(imag(v) - interval(im)))^2
            err2 += i == j ? e : 2e                # the mirrored entry has the same error
        end
    end
    return M, RationalBig(sup(sqrt(err2)))
end

"""
    verified_ring_spectrum(model, L; progress = false) -> Vector{BlockSpectrum}

Complete rigorous spectrum of the ring of L sites: every block of `ring_blocks`, enclosed by
`enclose_hermitian`.
"""
function verified_ring_spectrum(model::RingModel, L::Integer; progress = false)
    out = BlockSpectrum[]
    ring_blocks(model, L) do M, ε, mult, label
        seconds = @elapsed centers, radius = enclose_hermitian(M; matrix_error = ε)
        label["seconds"] = seconds
        push!(out, BlockSpectrum(centers, radius, mult, label))
        progress && println("L=$L charge=$(label["charge"]) m=$(label["momentum"]) ",
            "dim=$(size(M, 1)) radius=$(Float64(radius)) $(round(seconds; digits = 2))s")
    end
    return out
end

"""
    log_partition_enclosure(spectrum, β) -> Interval{Float64}

Enclosure of log Z(β) = log Σ e^{−βλ} over a verified spectrum, for exact β > 0.
"""
function log_partition_enclosure(spectrum::AbstractVector{BlockSpectrum}, β)
    b = to_interval(β)
    @assert inf(b) > 0
    shift = minimum(minimum(s.centers) for s in spectrum)
    total = interval(0.0)
    for s in spectrum
        part = sum(exp(-b * (interval(c) - interval(shift))) for c in s.centers)
        r = sup(to_interval(s.radius))
        total += s.multiplicity * part * exp(b * interval(-r, r))
    end
    return log(total) - b * interval(shift)
end

"JSON-ready summary of a verified spectrum (block labels, sizes, radii; no eigenvalues)."
function spectrum_summary(spectrum::AbstractVector{BlockSpectrum})
    Dict{String,Any}(
        "dimension" => sum(length(s.centers) * s.multiplicity for s in spectrum),
        "blocks" => length(spectrum),
        "largest_block" => maximum(length(s.centers) for s in spectrum),
        "max_radius" => string(maximum(s.radius for s in spectrum)),
        "max_radius_decimal" => Float64(maximum(s.radius for s in spectrum)),
        "seconds" => sum(s.label["seconds"] for s in spectrum),
    )
end

end
