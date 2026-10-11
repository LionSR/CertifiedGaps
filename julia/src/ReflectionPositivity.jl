"""
    CertifiedGaps.ReflectionPositivity

Reflection positivity of a bond in a chosen local basis (survey §1.2, Proposition 3.2). The
continuous-time transfer kernel is self-adjoint when −h = Σ_α C_α ⊗ conj(C_α) with a real
on-site term; equivalently the realignment R[−h]_{(ij),(kl)} = −h_{(ik),(jl)} is Hermitian and
positive semidefinite. Bond matrices act on `kron(site i, site j)`.
"""
module ReflectionPositivity

using LinearAlgebra

export realign, rp_test, rp_decomposition

"Realignment R[h]_{(i,j),(k,l)} = h_{(i,k),(j,l)} of a d²×d² bond matrix in kron order."
function realign(h::AbstractMatrix, d::Integer)
    # kron(A, B)[(i-1)d + k, (j-1)d + l] = A[i,j] B[k,l]; arrays are column-major, so
    # reshape(h, d, d, d, d)[k, i, l, j] = h[(i-1)d + k, (j-1)d + l].
    h4 = reshape(Matrix(h), d, d, d, d)
    R = Array{eltype(h)}(undef, d, d, d, d)          # R[(i,j),(k,l)] stored as R[j, i, l, k]
    for i in 1:d, j in 1:d, k in 1:d, l in 1:d
        R[j, i, l, k] = h4[k, i, l, j]
    end
    return reshape(R, d^2, d^2)
end

"""
    rp_test(h, d) -> (min_eig, antihermitian_norm)

Smallest eigenvalue of the Hermitian part of R[−h] and the norm of its anti-Hermitian part.
The bond is reflection positive in this basis iff the second is 0 and the first is ≥ 0.
"""
function rp_test(h::AbstractMatrix, d::Integer)
    R = realign(-complex(Matrix(h)), d)
    anti = norm((R - R') / 2)
    return eigmin(Hermitian((R + R') / 2)), anti
end

"Matrices C_α with −h = Σ_α C_α ⊗ conj(C_α), or `nothing` if the bond is not RP in this basis."
function rp_decomposition(h::AbstractMatrix, d::Integer; tol = 1e-10)
    R = realign(-complex(Matrix(h)), d)
    norm(R - R') > tol && return nothing
    w, v = eigen(Hermitian((R + R') / 2))
    minimum(w) < -tol && return nothing
    # R[A ⊗ conj(A)] = vec_row(A) vec_row(A)†, with vec_row the row-major vectorisation.
    [sqrt(λ) * permutedims(reshape(v[:, i], d, d)) for (i, λ) in enumerate(w) if λ > tol]
end

end
