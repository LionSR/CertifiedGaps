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
* or the polar factor U of X = UP: ‖X − U‖ ≤ ‖G‖/(1 + √(1 − ‖G‖)), and Weyl for
  U†MU = Λ + U†(MU − UΛ); the smaller of the two radii is kept;
* Weyl again for the distance between M and the true block.

Hence every true eigenvalue lies within one radius r of the sorted computed values. No
inequality rests on an unverified floating-point result. Scalar bounds are exact rationals.

`verified_eigen` (and its short form `enclose_hermitian`) is the enclosure step shared by every
model. A model with its own exact symmetry reduction (spin flip, reflection, leg exchange, …)
builds the Float64 center of each block and an exact bound on its distance to the true block
and passes both here; `ring_blocks` is the generic charge-and-momentum reduction.
"""
module VerifiedSpectra

using LinearAlgebra, IntervalArithmetic
using ..ExactArithmetic: RationalBig, parse_rational
using ..Results: write_result, read_result
using ..RingModels: RingModel, isreal_model, _reduce, _monomial
using ..Recurrence: to_interval
using ..ExactDiagonalization: PartialSpectrum

export BlockSpectrum, verified_eigen, enclose_hermitian, ring_blocks, verified_ring_spectrum,
    shifted_log_partition, log_partition_enclosure, ground_energy_enclosure,
    spectral_gap_enclosure, write_spectrum, read_spectrum, spectrum_summary, float_ring_spectrum

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
_sqrt_down(q::RationalBig) = RationalBig(inf(sqrt(interval(BigFloat, q))))

"""
    verified_eigen(M; matrix_error = 0) -> (; values, radius, residual, orthogonality, matrix_error)

Sorted eigenvalue approximations `values` of the exactly Hermitian Float64 matrix `M` (real or
complex) and an exact rational `radius` such that the sorted eigenvalues of every Hermitian A
with ‖A − M‖₂ ≤ `matrix_error` lie within `radius` of `values`, index by index. `residual` and
`orthogonality` are the exact upper bounds on ‖MX − XΛ‖₂ and ‖X†X − I‖₂ for the computed
eigenvectors X; `matrix_error` is returned as an exact rational.
"""
function verified_eigen(M::Matrix{T}; matrix_error = 0) where {T<:Union{Float64,ComplexF64}}
    ishermitian(M) || throw(ArgumentError("the block center is not exactly Hermitian"))
    N = size(M, 1)
    F = eigen(Hermitian(M))
    Λ, X = F.values, F.vectors
    @assert all(isfinite, Λ) && all(isfinite, X) && issorted(Λ)
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
    # Ostrowski: |λᵢ(M) − Λᵢ| ≤ (‖G‖ max|Λ| + ‖X‖ ‖R‖ + ‖G‖ |Λᵢ|)/(1 − ‖G‖), ‖X‖ ≤ 1 + ‖G‖.
    ostrowski = (2η * λmax + (1 + η) * r) / (1 - η)
    # Polar factor: the singular values of X lie in [√(1 − η), √(1 + η)], so ‖X − U‖ ≤ p. Then
    # U†MU − Λ = U†(MU − UΛ) is Hermitian with norm ≤ r + (‖M‖ + max|Λ|) p, and the first
    # bound gives ‖M‖ ≤ max|Λ| + `ostrowski`. This roughly halves the orthogonality term.
    p = η / (1 + _sqrt_down(1 - η))
    polar = r + (2λmax + ostrowski) * p
    return (; values = Λ, radius = min(ostrowski, polar) + RationalBig(matrix_error),
        residual = r, orthogonality = η, matrix_error = RationalBig(matrix_error))
end

"""
    enclose_hermitian(M; matrix_error = 0) -> (centers, radius)

The sorted `values` and the `radius` of `verified_eigen(M; matrix_error)`.
"""
function enclose_hermitian(M::Matrix{<:Union{Float64,ComplexF64}}; matrix_error = 0)
    e = verified_eigen(M; matrix_error)
    return e.values, e.radius
end

# State codes: digit of site x (0-based) is (code ÷ d^x) % d. The translation T moves the
# digit at site x to site x+1 (mod L).
_translate(code, d, top) = (code % top) * d + code ÷ top

"""
    ring_blocks(f, model, L; exact = true)

Exact symmetry decomposition of the ring of L sites into charge sectors, lattice momenta and the
parities of the model's on-site Z₂ `symmetries`. For each block it calls
`f(M, matrix_error, multiplicity, label)` with the Float64 center `M` of the interval block
matrix and an exact bound `matrix_error` on ‖A − M‖₂ for the true block A. Blocks are built
one at a time, so only one is in memory. For a real model the momenta ±k have the same
spectrum (A₋ₖ = conj Aₖ), so only 0 ≤ m ≤ L/2 is built, with multiplicity 2 for 0 < m < L/2.
For a model with a charge `flip`, the sectors q and −q have the same spectrum, so only one of
them is built (the one with the smaller |q|, else q > 0), with doubled multiplicity; a sector
with q ≡ −q is split by the parity of F^{⊗L} when the flip is a signed permutation commuting
with the symmetries. A model with a lattice `reflection` R also has equal spectra at ±k, so
only 0 ≤ m ≤ L/2 is built; for a real model with a reflection, the blocks at 0 < m < L/2 are
built in a real basis (see `_real_basis`), so they are real symmetric matrices.

The basis of a block is the symmetrized orbit states P_χ|r⟩ of the abelian group generated by
the translation T and the global parities G^{⊗L}, all acting on product states as signed
permutations. P_χ U_γ = χ(γ) P_χ, so if H r = Σ c_s |s⟩ with |s⟩ = σ U_γ |r_s⟩, then
⟨s̃|H|r̃⟩ = √(o_r/o_s) Σ c_s σ χ(γ), with o the orbit sizes; r contributes to the block of
character χ exactly when χ(δ) = σ for every element U_δ r = σ r of its stabilizer.

With `exact = false` the blocks are computed in plain Float64 arithmetic and `matrix_error` is
`nothing`: much faster, for floating-point scans only.
"""
function ring_blocks(f, model::RingModel, L::Integer; exact::Bool = true)
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
    reflection = model.reflection === nothing ? nothing : _monomial(model.reflection)
    paired = real_model || reflection !== nothing           # equal spectra at ±k
    momenta = paired ? (0:(L÷2)) : (0:(L-1))
    symmetries = [_monomial(G) for G in model.symmetries]
    flip = model.flip === nothing ? nothing : _monomial(model.flip)
    # The flip parity needs F^{⊗L} to be an involutive signed permutation in the abelian group.
    flip_parity = flip !== nothing && model.flip * model.flip == I &&
                  all(model.flip * G == G * model.flip for G in model.symmetries)
    total = 0
    for q in sort!(collect(keys(codes_by_charge)))
        partner = _reduce(model, -q)
        pair = model.flip === nothing || partner == q ? 1 : 2
        # Of q and −q, build the one with the smaller |q|, and q > 0 over −q.
        pair == 2 && (abs(partner), -partner) < (abs(q), -q) && continue
        maps = flip_parity && partner == q ? [symmetries; [flip]] : symmetries
        table = _sector_table(model, L, codes_by_charge[q], top, maps, reflection)
        for m in momenta, parities in Iterators.product(fill((1, -1), length(maps))...)
            built = _momentum_block(table, L, m, collect(Int, parities), real_model,
                _Arith{exact}())
            built === nothing && continue
            M, ε = built
            mult = pair * (paired && 0 < 2m < L ? 2 : 1)
            total += size(M, 1) * mult
            label = Dict{String,Any}("charge" => q, "momentum" => m, "dimension" => size(M, 1))
            isempty(symmetries) || (label["parities"] = collect(parities[1:length(symmetries)]))
            length(maps) > length(symmetries) && (label["flip_parity"] = last(parities))
            f(M, ε, mult, label)
        end
    end
    total == n || error("blocks cover $total of $n states")
    return nothing
end

# A global signed permutation ⊗ G (from `_monomial`) applied to a state code: (code′, sign)
# with G^{⊗L}|code⟩ = sign |code′⟩.
function _apply(code, (perm, sign), d, L, powers)
    out, σ, c = 0, 1, code
    for x in 1:L
        a = c % d + 1
        c ÷= d
        out += perm[a] * powers[x]
        σ *= sign[a]
    end
    return out, σ
end

# Orbit representatives of one sector under the group generated by T and the global signed
# permutations `maps` (bit j of g is the j-th map), and the exact action of H on them: for
# representative r, the coefficients c σ with H r = Σ c |s⟩ and |s⟩ = σ T^l G^g |r_s⟩, keyed
# by (index of r_s, l, g). Each representative carries its orbit size and the generators
# (l, g, σ) of its stabilizer, T^l G^g r = σ r, found as the cycles of the breadth-first search.
# With a `reflection` (from `_monomial`), `mirror[k]` is (k′, l, g, σ) with R r = σ T^l G^g r_k′
# for the lattice reflection R = P^{⊗L} ∘ (x ↦ L − 1 − x); otherwise `mirror` is `nothing`.
function _sector_table(model::RingModel, L, codes, top, maps, reflection = nothing)
    d = model.d
    powers = [d^x for x in 0:(L-1)]
    rep = Dict{Int,NTuple{4,Int}}()              # code => (representative index, l, g, σ)
    reps, sizes, stabilizers = Int[], Int[], Vector{Vector{NTuple{3,Int}}}()
    for code in codes                            # ascending, so each representative is minimal
        haskey(rep, code) && continue
        push!(reps, code)
        k = length(reps)
        rep[code] = (k, 0, 0, 1)
        queue, stabilizer, size = [code], Set{NTuple{3,Int}}(), 1
        while !isempty(queue)
            s = pop!(queue)
            _, l, g, σ = rep[s]
            for j in 0:length(maps)
                t, τ = j == 0 ? (_translate(s, d, top), 1) : _apply(s, maps[j], d, L, powers)
                l2, g2, σ2 = j == 0 ? mod(l + 1, L) : l, j == 0 ? g : g ⊻ (1 << (j - 1)), τ * σ
                if haskey(rep, t)
                    _, l3, g3, σ3 = rep[t]
                    push!(stabilizer, (mod(l3 - l2, L), g3 ⊻ g2, σ3 * σ2))
                else
                    rep[t] = (k, l2, g2, σ2)
                    push!(queue, t)
                    size += 1
                end
            end
        end
        push!(sizes, size)
        push!(stabilizers, collect(stabilizer))
    end
    T = isreal_model(model) ? Rational{Int} : Complex{Rational{Int}}
    convert_entry(x) = T === Rational{Int} ? Rational{Int}(real(x)) :
                       Complex{Rational{Int}}(Rational{Int}(real(x)), Rational{Int}(imag(x)))
    bond_cols = [[(i, convert_entry(model.bond[i, j])) for i in 1:d^2 if !iszero(model.bond[i, j])]
                 for j in 1:d^2]
    site_cols = [[(i, convert_entry(model.onsite[i, j])) for i in 1:d
                  if !iszero(model.onsite[i, j])]
                 for j in 1:d]
    columns = Vector{Dict{NTuple{3,Int},T}}(undef, length(reps))
    for (k, r) in enumerate(reps)
        column = Dict{NTuple{3,Int},T}()
        digit(x) = (r ÷ powers[x+1]) % d
        function add!(s, c)
            index, l, g, σ = rep[s]
            column[(index, l, g)] = get(column, (index, l, g), zero(T)) + σ * c
        end
        for x in 0:(L-1)
            y = (x + 1) % L
            a, b = digit(x), digit(y)
            for (i, c) in bond_cols[a*d+b+1]
                a2, b2 = divrem(i - 1, d)
                add!(r + (a2 - a) * powers[x+1] + (b2 - b) * powers[y+1], c)
            end
            for (i, c) in site_cols[a+1]
                add!(r + (i - 1 - a) * powers[x+1], c)
            end
        end
        columns[k] = column
    end
    mirror = reflection === nothing ? nothing : map(reps) do r
        (perm, sign), out, σ, c = reflection, 0, 1, r
        for x in 0:(L-1)                         # digit of site x goes to site L − 1 − x
            a = c % d + 1
            c ÷= d
            out += perm[a] * powers[L-x]
            σ *= sign[a]
        end
        k, l, g, τ = rep[out]
        (k, l, g, σ * τ)
    end
    return (; sizes, stabilizers, columns, mirror)
end

const _CI = Complex{Interval{Float64}}

# Arithmetic of the block entries: intervals (`_Arith{true}`, for the verified enclosure) or
# plain Float64 (`_Arith{false}`, for floating-point scans, which need no error bound).
struct _Arith{exact} end
_zero(::_Arith{true}) = complex(interval(0.0), interval(0.0))
_zero(::_Arith{false}) = 0.0im
_cispi(::_Arith{true}, x::Rational) =
    complex(cospi(interval(Float64, x)), sinpi(interval(Float64, x)))
_cispi(::_Arith{false}, x::Rational) = cispi(Float64(x))
_sqrt_ratio(::_Arith{true}, a, b) = sqrt(interval(Float64, a) / interval(Float64, b))
_sqrt_ratio(::_Arith{false}, a, b) = sqrt(a / b)
_number(::_Arith{true}, c) =
    complex(to_interval(Rational{BigInt}(real(c))), to_interval(Rational{BigInt}(imag(c))))
_number(::_Arith{false}, c) = ComplexF64(c)
_center(::_Arith{true}, x) = mid(x)
_center(::_Arith{false}, x) = x

# Product of the parities selected by the bits of g.
_parity(parities, g) = prod((parities[j] for j in eachindex(parities) if isodd(g >> (j - 1)));
    init = 1)

# Exact test of χ(l, g) = σ for the character χ(l, g) = e^{2πiml/L} · parity(g).
function _character_is(σ, L, m, parities, l, g)
    ε = σ * _parity(parities, g)                 # needs e^{2πiml/L} = ε ∈ {±1}
    return ε == 1 ? (m * l) % L == 0 : (2m * l) % L == 0 && (m * l) % L != 0
end

# Interval block at momentum k = 2πm/L and the given parities, then its center and an exact
# bound on ‖A − M‖_F:  ⟨s̃|H|r̃⟩ = √(o_r/o_s) Σ c σ e^{ikl} parity(g)  (see `ring_blocks`).
# With `_Arith{false}` the block is computed in Float64 and the bound is `nothing`.
function _momentum_block(table, L, m, parities, real_model, A = _Arith{true}())
    allowed = findall(eachindex(table.sizes)) do r
        all(_character_is(σ, L, m, parities, l, g) for (l, g, σ) in table.stabilizers[r])
    end
    isempty(allowed) && return nothing
    index = zeros(Int, length(table.sizes))
    index[allowed] = 1:length(allowed)
    N = length(allowed)
    exact = A isa _Arith{true}
    phases = [_cispi(A, 2 * ((m * l) % L) // L) for l in 0:(L-1)]
    ratios = Dict{Tuple{Int,Int},typeof(real(_zero(A)))}()
    coefficient = Dict{Any,typeof(_zero(A))}()
    enclose(c) = get!(() -> _number(A, c), coefficient, c)
    if real_model && table.mirror !== nothing && 0 < 2m < L
        return _real_block(table, L, m, parities, allowed, phases, ratios, enclose, A)
    end
    realblock = real_model && (2m) % L == 0
    M = realblock ? zeros(Float64, N, N) : zeros(ComplexF64, N, N)
    err2 = interval(0.0)
    for (j, r) in enumerate(allowed)
        entries = Dict{Int,typeof(_zero(A))}()
        for ((s, l, g), c) in table.columns[r]
            i = index[s]
            (i == 0 || i > j) && continue          # upper triangle, allowed states only
            ratio = get!(() -> _sqrt_ratio(A, table.sizes[r], table.sizes[s]), ratios,
                (table.sizes[r], table.sizes[s]))
            entries[i] = get(entries, i, _zero(A)) +
                         _parity(parities, g) * enclose(c) * phases[l+1] * ratio
        end
        for (i, v) in entries
            re = _center(A, real(v))
            im = (realblock || i == j) ? 0.0 : _center(A, imag(v))
            if realblock
                M[i, j] = re
                M[j, i] = re
            else
                M[i, j] = complex(re, im)
                M[j, i] = complex(re, -im)
            end
            exact || continue
            e = interval(mag(real(v) - interval(re)))^2 + interval(mag(imag(v) - interval(im)))^2
            err2 += i == j ? e : 2e                # the mirrored entry has the same error
        end
    end
    return M, exact ? RationalBig(sup(sqrt(err2))) : nothing
end

# Real basis of the block at 0 < k < π of a real model with a lattice reflection R. The
# antiunitary A = R ∘ conj commutes with H and fixes the momentum (R P_χ̄ R⁻¹ = P_χ, as
# R T R⁻¹ = T⁻¹ and R commutes with the parities), with A² = R² = 1, and acts on the orbit
# states as A|r̃⟩ = c|r̃′⟩, c = σ χ(l, g), for `mirror[r]` = (r′, l, g, σ). The A-fixed vectors
# (|r̃⟩ + c|r̃′⟩)/√2 and i(|r̃⟩ − c|r̃′⟩)/√2 of a pair r ≠ r′, and e^{iθ}|r̃⟩ with e^{2iθ} = c of a
# fixed r, form an orthonormal basis of the block in which ⟨u|H|v⟩ = conj⟨Au|H|Av⟩ is real.
# Each basis vector is a list of (orbit index, interval coefficient).
function _real_basis(table, L, m, parities, allowed, phases, A = _Arith{true}())
    inblock = falses(length(table.sizes))
    inblock[allowed] .= true
    h = _sqrt_ratio(A, 1, 2)
    basis = Vector{Vector{Tuple{Int,typeof(_zero(A))}}}()
    for r in allowed
        r2, l, g, σ = table.mirror[r]
        inblock[r2] && table.mirror[r2][1] == r || error("the reflection does not pair states")
        ε = σ * _parity(parities, g)
        c = ε * phases[l+1]
        if r2 == r
            # θ = π (ml mod L)/L, plus π/2 when ε = −1.
            push!(basis, [(r, _cispi(A, ((m * l) % L) // L + (ε == 1 ? 0 : 1 // 2)))])
        elseif r < r2
            push!(basis, [(r, complex(h, zero(h))), (r2, c * h)])
            push!(basis, [(r, complex(zero(h), h)), (r2, complex(zero(h), -h) * c)])
        end
    end
    return basis
end

# The block of `_momentum_block` in the real basis of `_real_basis`: a real symmetric center
# and an exact bound on ‖A − M‖_F. The true entries are real, so only the real part of each
# interval entry is kept; its imaginary part must contain zero (checked).
function _real_block(table, L, m, parities, allowed, phases, ratios, enclose, A)
    basis = _real_basis(table, L, m, parities, allowed, phases, A)
    T = typeof(_zero(A))
    members = [Tuple{Int,T}[] for _ in eachindex(table.sizes)]
    for (j, v) in enumerate(basis), (a, w) in v
        push!(members[a], (j, w))
    end
    exact = A isa _Arith{true}
    N = length(basis)
    M = zeros(Float64, N, N)
    err2 = interval(0.0)
    for (j, v) in enumerate(basis)
        image = Dict{Int,T}()                    # H v on the orbit states of the block
        for (r, w) in v, ((s, l, g), c) in table.columns[r]
            isempty(members[s]) && continue
            ratio = get!(() -> _sqrt_ratio(A, table.sizes[r], table.sizes[s]), ratios,
                (table.sizes[r], table.sizes[s]))
            image[s] = get(image, s, _zero(A)) +
                       w * (_parity(parities, g) * enclose(c) * phases[l+1] * ratio)
        end
        entries = Dict{Int,T}()
        for (s, y) in image, (i, w) in members[s]
            i <= j && (entries[i] = get(entries, i, _zero(A)) + conj(w) * y)
        end
        for (i, x) in entries
            re = _center(A, real(x))
            M[i, j] = re
            M[j, i] = re
            exact || continue
            in_interval(0, imag(x)) || error("a real-basis block entry is not real")
            e = interval(mag(real(x) - interval(re)))^2
            err2 += i == j ? e : 2e
        end
    end
    return M, exact ? RationalBig(sup(sqrt(err2))) : nothing
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

_interval(x::Interval) = x
_interval(x::Float64) = interval(x)
_interval(x) = to_interval(x)

"""
    shifted_log_partition(spectrum, β, shift) -> Interval{Float64}

Enclosure of log Σ e^{−β(λ − shift)} over a verified spectrum (multiplicities included), for β ≥ 0
exact or an `Interval` and an exact `shift`. Without the −β·shift term, two rings combined with
one exact shift difference keep a narrow enclosure over a whole β interval.
"""
function shifted_log_partition(spectrum::AbstractVector{BlockSpectrum}, β, shift)
    b, s0 = _interval(β), _interval(shift)
    inf(b) >= 0 || throw(ArgumentError("β must be nonnegative"))
    total = interval(0.0)
    for s in spectrum
        part = sum(exp(-b * (interval(c) - s0)) for c in s.centers)
        r = sup(to_interval(s.radius))
        total += s.multiplicity * part * exp(b * interval(-r, r))
    end
    return log(total)
end

"""
    log_partition_enclosure(spectrum, β) -> Interval{Float64}

Enclosure of log Z(β) = log Σ e^{−βλ} over a verified spectrum, for β ≥ 0 exact or an
`Interval`.
"""
function log_partition_enclosure(spectrum::AbstractVector{BlockSpectrum}, β)
    shift = minimum(minimum(s.centers) for s in spectrum)
    return shifted_log_partition(spectrum, β, shift) - _interval(β) * interval(shift)
end

# The two smallest eigenvalue balls of the ring, as (lower, upper) endpoint lists with
# multiplicity: every true eigenvalue lies in its ball, so the k-th smallest eigenvalue lies
# between the k-th smallest lower endpoint and the k-th smallest upper endpoint.
function _lowest_two(spectrum)
    lo, hi = RationalBig[], RationalBig[]
    for s in spectrum, c in s.centers[1:min(2, end)], _ in 1:min(2, s.multiplicity)
        push!(lo, RationalBig(c) - s.radius)
        push!(hi, RationalBig(c) + s.radius)
    end
    length(lo) >= 2 || throw(ArgumentError("the spectrum has fewer than two eigenvalues"))
    return partialsort!(lo, 1:2), partialsort!(hi, 1:2)
end

"""
    ground_energy_enclosure(spectrum) -> (lower, upper)

Exact rational bounds lower ≤ E₀ ≤ upper on the ground energy of a verified spectrum.
"""
function ground_energy_enclosure(spectrum::AbstractVector{BlockSpectrum})
    lo, hi = _lowest_two(spectrum)
    return lo[1], hi[1]
end

"""
    spectral_gap_enclosure(spectrum) -> (lower, upper)

Exact rational bounds lower ≤ E₁ − E₀ ≤ upper on the gap of a verified spectrum (E₁ is the second
eigenvalue counted with multiplicity, so a degenerate ground state has gap 0). `lower` is
clamped at 0; `lower > 0` proves a unique ground state.
"""
function spectral_gap_enclosure(spectrum::AbstractVector{BlockSpectrum})
    lo, hi = _lowest_two(spectrum)
    return max(lo[2] - hi[1], zero(RationalBig)), hi[2] - lo[1]
end

const SPECTRUM_FORMAT = "CertifiedGaps verified spectrum v1"

"""
    write_spectrum(path, spectrum) -> path

Store a verified spectrum losslessly (`.json` or `.json.gz` via `Results.write_result`): exact
radii, multiplicities, labels, and every center as its shortest round-trip decimal. The output
is deterministic, so its SHA-256 can pin the spectrum in a later certificate. The file is read
back and compared before returning.
"""
function write_spectrum(path, spectrum::AbstractVector{BlockSpectrum})
    blocks = [Dict{String,Any}("centers" => string.(s.centers), "radius" => string(s.radius),
        "multiplicity" => s.multiplicity, "label" => s.label) for s in spectrum]
    write_result(path, Dict{String,Any}("format" => SPECTRUM_FORMAT, "blocks" => blocks);
        compact = true)
    back = read_spectrum(path)
    all(a.centers == b.centers && a.radius == b.radius && a.multiplicity == b.multiplicity
        for (a, b) in zip(back, spectrum)) && length(back) == length(spectrum) ||
        error("$path does not read back as the spectrum written")
    return path
end

"Inverse of `write_spectrum`."
function read_spectrum(path)
    data = read_result(path; exact = false)
    data["format"] == SPECTRUM_FORMAT || error("$path is not a verified spectrum file")
    return [BlockSpectrum(parse.(Float64, b["centers"]), parse_rational(b["radius"]),
        Int(b["multiplicity"]), Dict{String,Any}(String(k) => copy(v) for (k, v) in b["label"]))
            for b in data["blocks"]]
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

"""
    float_ring_spectrum(model, L) -> PartialSpectrum

All d^L levels of the ring of L sites, sorted, from unverified LAPACK eigenvalues of the
`ring_blocks` blocks with their exact multiplicities. Floating point only: for seed scans and
diagnostics, never a certificate (`verified_ring_spectrum` is the verified counterpart).
"""
function float_ring_spectrum(model::RingModel, L::Integer)
    E = Float64[]
    ring_blocks(model, L; exact = false) do M, _, mult, _
        e = eigvals(Hermitian(M))
        for _ in 1:mult
            append!(E, e)
        end
    end
    return PartialSpectrum(E, model.d^L)
end

end
