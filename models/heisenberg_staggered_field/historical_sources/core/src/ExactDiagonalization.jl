"""
    ThermalMPOGap.ExactDiagonalization

Exact diagonalisation of spin systems with nearest-neighbour or graph interactions, split into
sectors of a conserved charge, and partition-function enclosures from the resulting spectra
(survey Proposition 4.12).

Conventions: basis states are base-d integers with **site 1 the most significant digit**; a
two-site term `h` acts on `kron(site i, site j)`. Spectra do not depend on these choices.
"""
module ExactDiagonalization

using KrylovKit, LinearAlgebra, SparseArrays

export Chain, ring_bonds, open_bonds, graph_bonds, sectors, sector_states,
    sector_hamiltonian, PartialSpectrum, spectrum, logZ_bounds, log_partition

"""
    Chain(d, bonds; onsite = nothing, charge = nothing)

A system of local dimension `d`. `bonds(L)` returns the two-site terms of the system of `L`
sites as `(i, j, h)` triples, `h` a d²×d² matrix on `kron(site i, site j)`. `onsite` is a d×d
term on every site. `charge` is an integer per local state; when given, `H` must conserve
Σ_j charge[s_j], and spectra are computed sector by sector.
"""
struct Chain{B,O,Q}
    d::Int
    bonds::B
    onsite::O
    charge::Q
end
Chain(d, bonds; onsite = nothing, charge = nothing) = Chain(d, bonds, onsite, charge)

"Periodic ring: bond j joins sites j and j+1 mod L (for L = 2 both bonds join sites 1, 2)."
ring_bonds(h) = L -> [(j, mod1(j + 1, L), h) for j in 1:L]

"Open chain: bonds (j, j+1) for j < L."
open_bonds(h) = L -> [(j, j + 1, h) for j in 1:(L-1)]

"The same term `h` on every pair of a fixed list (an arbitrary graph)."
graph_bonds(pairs, h) = _ -> [(i, j, h) for (i, j) in pairs]

"Digit (0-based local state) of site j in basis state s."
digit(s, j, L, d) = (s ÷ d^(L - j)) % d

"""
    sectors(c, L) -> Dict(charge => sorted states)

All basis states of `L` sites grouped by total charge (one group, keyed `nothing`, if the
chain has no charge).
"""
function sectors(c::Chain, L::Integer)
    N = c.d^L
    c.charge === nothing && return Dict{Any,Vector{Int}}(nothing => collect(0:(N-1)))
    groups = Dict{Int,Vector{Int}}()
    for s in 0:(N-1)
        q = sum(c.charge[digit(s, j, L, c.d)+1] for j in 1:L)
        push!(get!(groups, q, Int[]), s)
    end
    return groups
end

"Sorted basis states of `L` sites with total charge `q` (all states if `q === nothing`)."
sector_states(c::Chain, L::Integer, q = nothing) =
    q === nothing || c.charge === nothing ? collect(0:(c.d^L-1)) : get(sectors(c, L), q, Int[])

"""
    sector_hamiltonian(c, L, states) -> SparseMatrixCSC

The Hamiltonian of `L` sites restricted to the span of `states` (which must be closed under
`H`, e.g. a charge sector).
"""
function sector_hamiltonian(c::Chain, L::Integer, states::AbstractVector{<:Integer})
    d = c.d
    index = Dict(s => i for (i, s) in enumerate(states))
    place = [d^(L - j) for j in 1:L]
    I, J, V = Int[], Int[], ComplexF64[]
    terms = _local_terms(c, L)
    for (col, s) in enumerate(states), term in terms
        old = [digit(s, j, L, d) for j in term.sites]
        oldidx = _local_index(old, d)
        for newidx in axes(term.mat, 1)
            amp = term.mat[newidx, oldidx]
            abs(amp) > 1e-14 || continue
            new = _local_digits(newidx, length(old), d)
            t = s + sum((new[k] - old[k]) * place[term.sites[k]] for k in eachindex(old))
            row = get(index, t, 0)
            row == 0 && error("term leaves the sector: charge not conserved")
            push!(I, row)
            push!(J, col)
            push!(V, amp)
        end
    end
    H = sparse(I, J, V, length(states), length(states))
    return maximum(abs, imag(H); init = 0.0) < 1e-14 ? real(H) : H
end

const LocalTerm = @NamedTuple{sites::Vector{Int}, mat::Matrix{ComplexF64}}

# Two-site terms first, then on-site terms. Computed once per Hamiltonian.
function _local_terms(c::Chain, L)
    terms = LocalTerm[(sites = [i, j], mat = h) for (i, j, h) in c.bonds(L)]
    c.onsite === nothing || append!(terms, LocalTerm[(sites = [j], mat = c.onsite) for j in 1:L])
    return terms
end

_local_index(digits, d) = foldl((acc, x) -> acc * d + x, digits; init = 0) + 1
_local_digits(idx, n, d) = [((idx - 1) ÷ d^(n - k)) % d for k in 1:n]

"""
Lowest eigenvalues `E` (sorted, with multiplicity) of an operator of dimension `N`; `exact`
when `E` is the whole spectrum. Every uncomputed eigenvalue is at least `E[end]`.
"""
struct PartialSpectrum
    E::Vector{Float64}
    N::Int
    exact::Bool
end
PartialSpectrum(E, N) = PartialSpectrum(sort(E), N, length(E) == N)

"""
    spectrum(c, L; spin_flip = false, dense_max = 3000, k = 60, on_sector = nothing)

The spectrum of `L` sites, merged over charge sectors. Sectors up to `dense_max` states are
diagonalised densely (exact); larger ones by Lanczos (lowest `k` levels). Lanczos can miss
degenerate partners, so the result is then not a certificate. When a sector is truncated, all
levels are cut at the lowest top-kept level across sectors, which keeps the partition-function
enclosure valid.

With `spin_flip = true` the chain must be symmetric under charge → −charge; only sectors with
charge ≥ 0 are diagonalised and those with charge > 0 are counted twice.

`on_sector(row)` is called after each sector with
`row = (; charge, dimension, multiplicity, levels, seconds)`.
"""
function spectrum(c::Chain, L::Integer; spin_flip = false, dense_max = 3000, k = 60,
    on_sector = nothing)
    groups = sectors(c, L)
    spin_flip && _check_flip_symmetric(c)
    charges = sort!(collect(keys(groups)); by = q -> q === nothing ? 0 : q)
    spin_flip && filter!(q -> q >= 0, charges)
    levels, N, cut = Float64[], 0, Inf
    for q in charges
        start = time()
        E, complete = _low_levels(sector_hamiltonian(c, L, groups[q]), dense_max, k)
        complete || (cut = min(cut, E[end]))
        multiplicity = spin_flip && q > 0 ? 2 : 1
        for _ in 1:multiplicity
            append!(levels, E)
        end
        N += multiplicity * length(groups[q])
        on_sector === nothing || on_sector((; charge = q, dimension = length(groups[q]),
            multiplicity, levels = E, seconds = time() - start))
    end
    isfinite(cut) && filter!(≤(cut), levels)
    return PartialSpectrum(sort!(levels), N, isinf(cut))
end

function _check_flip_symmetric(c::Chain)
    c.charge === nothing && throw(ArgumentError("spin_flip needs a charge vector"))
    c.charge == -reverse(c.charge) ||
        throw(ArgumentError("spin_flip needs charge → −charge symmetry"))
end

function _low_levels(H, dense_max, k)
    n = size(H, 1)
    n <= dense_max && return eigvals(Hermitian(Matrix(H))), true
    vals, _, info = eigsolve(H, rand(eltype(H), n), k, :SR;
        ishermitian = true, krylovdim = max(2k, 100), tol = 1e-11, maxiter = 500)
    info.converged >= k || @warn "Lanczos converged $(info.converged) of $k levels"
    return sort(real(vals[1:min(k, info.converged)])), false
end

"""
    log_partition(energies, β)

log Σ exp(−β E) over a complete list of levels, shifted by the minimum for stability.
"""
function log_partition(energies, beta)
    beta >= 0 || throw(ArgumentError("beta must be nonnegative"))
    e0 = minimum(energies)
    return -beta * e0 + log(sum(exp(-beta * (e - e0)) for e in energies))
end

"""
    logZ_bounds(spec, β) -> (lo, hi)

Enclosure of log Z(β): the kept levels give the lower end, and the N − k uncomputed levels,
each at least E[end], the upper end. Exact when the spectrum is complete.
"""
function logZ_bounds(spec::PartialSpectrum, β::Real)
    E₀ = spec.E[1]
    lo = _logsumexp(-β .* (spec.E .- E₀))
    hi = spec.exact ? lo :
         _logaddexp(lo, log(spec.N - length(spec.E)) - β * (spec.E[end] - E₀))
    return lo - β * E₀, hi - β * E₀
end

_logaddexp(a, b) = max(a, b) + log1p(exp(-abs(a - b)))
_logsumexp(x) = (m = maximum(x); m + log(sum(exp.(x .- m))))

end
