"Square-lattice Ising tori represented as rings of complete transverse rows."
module SquareIsing

using LinearAlgebra
using CertifiedGaps.RingModels: RingModel, rp_certificate, verify_rp

export row_model, row_operators, row_rp

"Exact Pauli X and Z acting at each position of a row of width m."
function row_operators(m::Integer)
    m in (2, 4) || throw(ArgumentError("the finite seeds use widths two and four"))
    d = 2^m
    xs = [zeros(Int, d, d) for _ in 1:m]
    zs = [zeros(Int, d, d) for _ in 1:m]
    for s in 0:d-1, y in 0:m-1
        xs[y + 1][xor(s, 1 << y) + 1, s + 1] = 1
        zs[y + 1][s + 1, s + 1] = iszero(s & (1 << y)) ? 1 : -1
    end
    return xs, zs
end

function exchange(m, a, b)
    d = 2^m
    p = zeros(Int, d, d)
    for s in 0:d-1
        t = ((s >> a) & 1) == ((s >> b) & 1) ? s : xor(s, (1 << a) | (1 << b))
        p[t + 1, s + 1] = 1
    end
    return p
end

"""
    row_model(m, g)

The Hadamard-rotated Hamiltonian −Σ XX − g Σ Z, as a ring of rows of width m.
The transverse periodic sum retains both directed occurrences when m = 2.
Opposite-vertex exchanges on the four-cycle commute and reduce the spectral blocks.
"""
function row_model(m::Integer, g)
    xs, zs = row_operators(m)
    d = 2^m
    bond = -sum(kron(x, x) for x in xs)
    onsite = -sum(xs[y] * xs[mod1(y + 1, m)] for y in 1:m) - g * sum(zs)
    syms = m == 2 ? [exchange(m, 0, 1)] : [exchange(m, 0, 2), exchange(m, 1, 3)]
    return RingModel(d, bond, onsite; charge = [count_ones(s) % 2 for s in 0:d-1],
        modulus = 2, symmetries = syms, reflection = Matrix{Int}(I, d, d),
        name = "square TFIM, transverse period $m, g = $g, J = 1")
end

"Exact row reflection-positive factorization, verified by the shared core."
function row_rp(m::Integer, g)
    model = row_model(m, g)
    xs, _ = row_operators(m)
    cert = rp_certificate(model; factors = xs)
    verify_rp(model, cert) || error("row reflection positivity failed")
    return cert
end

end
