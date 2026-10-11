module ThermalBounds

using LinearAlgebra
import Nemo

"Exact two-spin scalar Hamiltonian in a rational Marshall gauge."
function rational_bond()
    plus, minus = zeros(Int, 5, 5), zeros(Int, 5, 5)
    for r in 0:3
        plus[r+2, r+1] = 4-r
        minus[r+1, r+2] = r+1
    end
    z = Matrix(Diagonal(-2:2))
    h = Nemo.matrix(Nemo.QQ, 2kron(z, z)-kron(plus, minus)-kron(minus, plus))/2
    eye = Nemo.identity_matrix(Nemo.QQ, 25)
    eigenvalues = [-6, -5, -3, 0, 4]
    projectors = [prod((h-f*eye)/(e-f) for f in eigenvalues if f != e) for e in eigenvalues]
    @assert sum(projectors) == eye
    @assert sum(e*p for (e, p) in zip(eigenvalues, projectors)) == h
    for (i, p) in enumerate(projectors), (j, q) in enumerate(projectors)
        @assert p*q == (i == j ? p : zero(p))
    end
    return (; h, eigenvalues, projectors)
end

function gate_matrix(bond, delta, RR)
    gate =
        sum(exp(-delta*e)*Nemo.matrix(RR, p) for (e, p) in zip(bond.eigenvalues, bond.projectors))
    # Conservation makes these zeros exact, avoiding interval cancellation.
    for a in 0:24, b in 0:24
        div(a, 5)+mod(a, 5) == div(b, 5)+mod(b, 5) || (gate[a+1, b+1] = 0)
    end
    return gate
end

function layer_matrix(n, magnetization, gate, offset, RR)
    powers = [5^(j-1) for j in 1:n]
    states =
        [code for code in 0:(5^n-1) if sum(mod(div(code, p), 5)-2 for p in powers) == magnetization]
    layer = Nemo.zero_matrix(RR, length(states), length(states))
    pairs = [(j, mod1(j+1, n)) for j in offset:2:n]
    for (i, row) in enumerate(states), (j, col) in enumerate(states)
        value = RR(1)
        for (a, b) in pairs
            r = 5mod(div(row, powers[a]), 5)+mod(div(row, powers[b]), 5)+1
            c = 5mod(div(col, powers[a]), 5)+mod(div(col, powers[b]), 5)+1
            value *= gate[r, c]
        end
        layer[i, j] = value
    end
    return layer
end

"""
Enclose the Trotter trace upper bound using rigorous Arb arithmetic.

For even periodic n, H=A+B consists of two disjoint bond layers. ALT gives
Z(beta) <= Tr[(exp(-beta*A/m) exp(-beta*B/m))^m]. This needs no numerical
estimate of the Trotter error. The returned ball encloses the RIGHT-hand
side, and is not a two-sided enclosure of the true partition function.

This full-sector implementation is a reference verifier, limited to n<=6.
It has no MPO compression; scalable certified contraction remains necessary.
"""
function trotter_upper(n, beta, m; shift = 0, precision = 192)
    iseven(n) && 4 <= n <= 6 || throw(ArgumentError("reference lengths are 4 and 6"))
    m isa Integer && m >= 1 && beta > 0 ||
        throw(ArgumentError("positive beta and integer step count required"))
    beta isa Union{Integer,Rational} || throw(ArgumentError("beta must be exact"))
    shift isa Union{Integer,Rational} || throw(ArgumentError("shift must be exact"))
    RR = Nemo.ArbField(precision)
    bond = rational_bond()
    gate = gate_matrix(bond, RR(beta)/m, RR)
    total = RR(0)
    sectors = Any[]
    for mag in 0:2n
        a = layer_matrix(n, mag, gate, 1, RR)
        b = layer_matrix(n, mag, gate, 2, RR)
        value = Nemo.tr((a*b)^m)
        total += (mag == 0 ? 1 : 2)*value
        push!(sectors, (; magnetization = mag, dimension = size(a, 1), trace = value))
    end
    prefactor = exp(-RR(beta)*RR(shift)*n)
    return (;
        trotter_trace = prefactor*total,
        log_trotter_trace = log(total)-RR(beta)*RR(shift)*n,
        sectors,
    )
end

end
