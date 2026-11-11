# Sparse symmetry sectors of the physical-frame ladder ring of L rungs (N = 2L spins), floating
# point only (screening, not certified). Included by scripts/ritz_scan.jl.
#
# Spin (j, leg) is bit 2j + leg of a state code (leg a = 0, leg b = 1), Z = (−1)^bit. The
# abelian group G = ⟨T⟩ × ⟨S⟩ × ⟨F⟩ of order 4L (rung translation T, leg exchange S, global
# spin flip F = ∏ X) acts on codes as permutations, so the sector of the character
# χ(T^l S^a F^b) = e^{2πiml/L} s^a f^b is spanned by the normalized projections of the orbit
# representatives r whose stabilizer lies in ker χ, with
#     ⟨r̃′|H|r̃⟩ = √(o_r / o_r′) Σ c χ(γ)    for  H r = Σ c U_γ r′  (o = orbit size).
# The lattice reflection (not used here) makes momenta ±k isospectral, so only 0 ≤ m ≤ L/2
# is built, with multiplicity 2 for 0 < m < L/2.
using SparseArrays, LinearAlgebra

struct LadderGroup
    L::Int
    N::Int
    mask::UInt64
    even::UInt64
end
LadderGroup(L) = LadderGroup(L, 2L, (UInt64(1) << 2L) - 1,
    reduce(|, (UInt64(1) << 2j for j in 0:(L-1))))

order(G::LadderGroup) = 4G.L
# Group element e ∈ 0:4L−1 encodes (l, a, b) = (e ÷ 4, (e ÷ 2) % 2, e % 2).
@inline function act(G::LadderGroup, x::UInt64, e::Int)
    l, a, b = e >> 2, (e >> 1) & 1, e & 1
    if l > 0
        s = 2l
        x = ((x << s) | (x >> (G.N - s))) & G.mask
    end
    a == 1 && (x = ((x & G.even) << 1) | ((x >> 1) & G.even))
    b == 1 && (x = x ⊻ G.mask)
    return x
end

"(representative, e) with act(G, x, e) = representative, the minimal code of the orbit of x."
@inline function representative(G::LadderGroup, x::UInt64)
    best, be = x, 0
    for e in 1:(order(G)-1)
        y = act(G, x, e)
        y < best && ((best, be) = (y, e))
    end
    return best, be
end

"""
    orbit_table(G) -> (; reps, stab, size)

All orbit representatives (sorted), the stabilizer of each as a bit mask over group elements,
and the orbit sizes.
"""
function orbit_table(G::LadderGroup)
    reps, stab = UInt64[], UInt64[]
    n = order(G)
    for x in UInt64(0):G.mask
        minimal, mask = true, UInt64(1)
        for e in 1:(n-1)
            y = act(G, x, e)
            if y < x
                minimal = false
                break
            end
            y == x && (mask |= UInt64(1) << e)
        end
        minimal || continue
        push!(reps, x)
        push!(stab, mask)
    end
    return (; reps, stab, size = [n ÷ count_ones(s) for s in stab])
end

"Diagonal energy −Σ Z Z (legs and rungs) of a code."
@inline function zz_energy(G::LadderGroup, x::UInt64)
    rot = ((x << 2) | (x >> (G.N - 2))) & G.mask
    leg_anti = count_ones(x ⊻ rot)                       # anti-aligned leg bonds
    rung_anti = count_ones((x ⊻ (x >> 1)) & G.even)      # anti-aligned rungs
    return -(2G.L - 2leg_anti) - (G.L - 2rung_anti)
end

"""
    hopping_table(G, table) -> (target, element)

For every representative r and spin i, the representative index of the flipped state r ⊻ 2^i
and the group element e with act(e) mapping that state to its representative.
"""
function hopping_table(G::LadderGroup, table)
    nr = length(table.reps)
    target = Matrix{Int32}(undef, G.N, nr)
    element = Matrix{UInt8}(undef, G.N, nr)
    for k in 1:nr, i in 0:(G.N-1)
        r, e = representative(G, table.reps[k] ⊻ (UInt64(1) << i))
        target[i+1, k] = searchsortedfirst(table.reps, r)
        element[i+1, k] = e
    end
    return target, element
end

# χ(γ) for the element e and the character (m, s, f).
function character(G::LadderGroup, e, m, s, f)
    l, a, b = e >> 2, (e >> 1) & 1, e & 1
    return cispi(2m * l / G.L) * (a == 1 ? s : 1) * (b == 1 ? f : 1)
end

compatible(G, stab, m, s, f) = all(e -> !isodd(stab >> e) ||
    isapprox(character(G, e, m, s, f), 1; atol = 1e-12), 0:(order(G)-1))

"""
    sector_matrix(G, table, hops, g, m, s, f) -> SparseMatrixCSC

The sector block of the character (m, s, f) of H = −Σ ZZ − g Σ X, real for m ∈ {0, L/2}.
"""
function sector_matrix(G::LadderGroup, table, hops, g, m, s, f)
    target, element = hops
    nr = length(table.reps)
    idx = zeros(Int32, nr)
    dim = 0
    for k in 1:nr
        compatible(G, table.stab[k], m, s, f) && (idx[k] = (dim += 1))
    end
    real_sector = m == 0 || 2m == G.L
    T = real_sector ? Float64 : ComplexF64
    I, J, V = Int32[], Int32[], T[]
    sizehint!(I, dim * (G.N + 1)); sizehint!(J, dim * (G.N + 1)); sizehint!(V, dim * (G.N + 1))
    for k in 1:nr
        col = idx[k]
        col == 0 && continue
        push!(I, col); push!(J, col); push!(V, zz_energy(G, table.reps[k]))
        for i in 1:G.N
            t = target[i, k]
            row = idx[t]
            row == 0 && continue
            # H r ∋ −g |x⟩ with x = U_e⁻¹ r′, so χ(e⁻¹) = conj χ(e).
            c = -g * sqrt(table.size[k] / table.size[t]) *
                conj(character(G, element[i, k], m, s, f))
            push!(I, row); push!(J, col); push!(V, real_sector ? real(c) : c)
        end
    end
    return sparse(I, J, V, dim, dim)
end

"The characters (m, s, f) with 0 ≤ m ≤ L/2 and their multiplicities (2 for 0 < m < L/2)."
sectors(L) = [(m, s, f, (0 < 2m < L) ? 2 : 1) for m in 0:(L÷2) for s in (1, -1) for f in (1, -1)]
