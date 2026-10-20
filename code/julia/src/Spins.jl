"""
    CertifiedGaps.Spins

Local spin operators in the three bases the projects use:

* `spin_matrices(s)`: the standard basis m = s, s−1, …, −s (descending), normalised;
* `integer_ladders(j)`: the basis m = −j, …, j (ascending) in the unnormalised binomial gauge,
  with integer entries, for exact arithmetic;
* `cartesian_generators(s)`: a real basis in which G_α = iS^α are real antisymmetric, so that
  −S·S = Σ_α G_α ⊗ G_α is reflection positive (s = 1, 2);
* exact matrices for RP identities: `exact_pauli`, `exact_spin_half` (any rational element
  type), `integer_cartesian_generators` and `integer_spin_one` (spin 1, integer entries).
"""
module Spins

using LinearAlgebra

export spin_matrices, integer_ladders, binomial_weights, cartesian_generators, PAULI,
    exact_pauli, exact_spin_half, integer_cartesian_generators, integer_spin_one

"""
    spin_matrices(s) -> (; sz, sp, sm, sx, sy)

Spin-s operators in the basis m = s, s−1, …, −s. All are `Matrix{Float64}` except `sy`, which
is complex. Entry for entry equal to spin2's and aklt's `spin_matrices()`.
"""
function spin_matrices(s::Real)
    d = round(Int, 2s + 1)
    m = [s - (j - 1) for j in 1:d]
    sz = Matrix(Diagonal(Float64.(m)))
    sp = zeros(d, d)
    for j in 2:d
        sp[j-1, j] = sqrt(s * (s + 1) - m[j] * (m[j] + 1))
    end
    sm = copy(transpose(sp))
    return (; sz, sp, sm, sx = (sp + sm) / 2, sy = (sp - sm) / (2im))
end

"""
    integer_ladders(j) -> (plus, minus)

Raising and lowering operators for spin j in the basis m = −j, …, j, in the binomial gauge
where both are integer matrices: `plus[r+2, r+1] = 2j − r`, `minus[r+1, r+2] = r + 1`.
They satisfy `plus * Diagonal(w) == Diagonal(w) * transpose(minus)` with `w = binomial_weights(j)`.
"""
function integer_ladders(j)
    n = 2j + 1
    plus, minus = zeros(BigInt, n, n), zeros(BigInt, n, n)
    for r in 0:(2j-1)
        plus[r+2, r+1] = 2j - r
        minus[r+1, r+2] = r + 1
    end
    return plus, minus
end

"Metric of the binomial gauge: `binomial(2j, j + m)` for m = −j, …, j."
binomial_weights(j) = [binomial(Int(2j), k) for k in 0:Int(2j)]

"""
    cartesian_generators(s) -> [G_x, G_y, G_z]

Real antisymmetric G_α = iS^α for s = 1 (vectors, `(G_a)_{bc} = ε_{abc}`) and s = 2 (symmetric
traceless 3×3 matrices, Q ↦ J_α Q − Q J_α in an orthonormal basis). In this basis
−S·S = Σ_α G_α ⊗ G_α, which is reflection positive.
"""
function cartesian_generators(s::Integer)
    s == 1 && return [_rotation_generator(a) for a in 1:3]
    s == 2 && return _quadrupole_generators()
    throw(ArgumentError("Cartesian generators are implemented for s = 1, 2"))
end

# Real antisymmetric generators of rotations of ℝ³: (J_a)_{bc} = ε_{abc}.
_rotation_generator(a) = Float64.(integer_cartesian_generators()[a])

_levi_civita(a, b, c) = (a, b, c) in ((1, 2, 3), (2, 3, 1), (3, 1, 2)) ? 1 :
                        (a, b, c) in ((1, 3, 2), (3, 2, 1), (2, 1, 3)) ? -1 : 0

"""
    integer_cartesian_generators() -> [G_x, G_y, G_z]

Spin-1 Cartesian generators with integer entries, `(G_a)_{bc} = ε_{abc}` (`Matrix{Int}`): real
antisymmetric, `S^a = −i G_a` satisfies [S^x, S^y] = iS^z and Σ_a (S^a)² = 2, and
−S·S = Σ_a G_a ⊗ G_a. Equal to `cartesian_generators(1)` entry for entry.
"""
integer_cartesian_generators() = [[_levi_civita(a, b, c) for b in 1:3, c in 1:3] for a in 1:3]

# Moved verbatim from spin2 `model.jl` (`real_generators`). Norm 2, not the spin-1 norm 1.
function _quadrupole_generators()
    basis = [
        Float64[1 0 0; 0 -1 0; 0 0 0] / sqrt(2),
        Float64[1 0 0; 0 1 0; 0 0 -2] / sqrt(6),
        Float64[0 1 0; 1 0 0; 0 0 0] / sqrt(2),
        Float64[0 0 1; 0 0 0; 1 0 0] / sqrt(2),
        Float64[0 0 0; 0 0 1; 0 1 0] / sqrt(2),
    ]
    generators = [
        Float64[0 0 0; 0 0 1; 0 -1 0],
        Float64[0 0 -1; 0 0 0; 1 0 0],
        Float64[0 1 0; -1 0 0; 0 0 0],
    ]
    return [[dot(a, j * b - b * j) for a in basis, b in basis] for j in generators]
end

"""
    exact_pauli(T = Rational{BigInt}) -> (; x, y, z, iy, id)

Exact Pauli matrices with element type `T` (any exact real type): `x`, `z`, `id` and the real
antisymmetric `iy = iσʸ = [0 1; −1 0]` are `Matrix{T}`, `y = σʸ` is `Matrix{Complex{T}}`.
"""
exact_pauli(::Type{T} = Rational{BigInt}) where {T<:Real} = (x = T[0 1; 1 0],
    y = Complex{T}[0 -im; im 0], z = T[1 0; 0 -1], iy = T[0 1; -1 0], id = T[1 0; 0 1])

"""
    exact_spin_half(T = Rational{BigInt}) -> (; x, y, z, iy, id)

Exact spin-½ operators S = σ/2 in the basis (↑, ↓): `x`, `z`, the real antisymmetric
`iy = iSʸ = [0 1; −1 0]/2` and `id` are `Matrix{T}`, `y = Sʸ` is `Matrix{Complex{T}}`. `T` must
represent halves exactly (a rational type). With real factors, S·S = x⊗x − iy⊗iy + z⊗z.
"""
function exact_spin_half(::Type{T} = Rational{BigInt}) where {T<:Real}
    σ = exact_pauli(T)
    return (x = σ.x / 2, y = σ.y / 2, z = σ.z / 2, iy = σ.iy / 2, id = σ.id)
end

"Pauli matrices and the identity."
const PAULI = (x = ComplexF64[0 1; 1 0], y = ComplexF64[0 -im; im 0],
    z = ComplexF64[1 0; 0 -1], i = ComplexF64[1 0; 0 1])

"""
    integer_spin_one() -> (; xp, ap, sz, flip, id)

Spin 1 in the Sᶻ basis m = 1, 0, −1 with integer entries (`Matrix{Int}`): `xp = √2 Sˣ`, the
real antisymmetric `ap = √2 iSʸ`, `sz = Sᶻ`, and the spin flip `flip` (m ↦ −m, which keeps Sˣ
and reverses Sʸ and Sᶻ). Then S·S = xp⊗xp/2 − ap⊗ap/2 + sz⊗sz exactly.
"""
integer_spin_one() = (xp = [0 1 0; 1 0 1; 0 1 0], ap = [0 1 0; -1 0 1; 0 -1 0],
    sz = [1 0 0; 0 0 0; 0 0 -1], flip = [0 0 1; 0 1 0; 1 0 0], id = [1 0 0; 0 1 0; 0 0 1])

end
