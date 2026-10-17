# Floating-point spectra for the seed scan (`scripts/scan.jl`), never an input to a certificate:
# the physical S^z basis, exact integer Hamiltonian blocks, and the real dihedral × spin-flip
# reduction (cell translation, site reversal, global flip), which keeps the 12-spin ring at
# blocks of dimension at most 11 596. Part of module `DimerizedCertificate`.
using SparseArrays
using CertifiedGaps.ExactDiagonalization: PartialSpectrum

"Exact spin-1 S^z-basis Hamiltonian of a block of product states, scaled to an integer matrix."
function integer_hamiltonian(cells, delta, states)
    L = 2cells; den = denominator(delta); num = numerator(delta)
    index = Dict(s => i for (i, s) in enumerate(states)); place = 3 .^ (0:L-1)
    rows = Int[]; cols = Int[]; vals = Int[]
    for (col, s) in enumerate(states)
        diag = 0
        for j in 1:L
            k = mod1(j + 1, L); a = (s ÷ place[j]) % 3; b = (s ÷ place[k]) % 3
            w = den + (iseven(j) ? num : -num)
            diag += w * (a - 1) * (b - 1)
            # Spin-1 exchange: the two √2 ladder factors cancel the 1/2 of (S⁺S⁻ + S⁻S⁺)/2.
            for e in (-1, 1)
                if 0 <= a + e <= 2 && 0 <= b - e <= 2
                    push!(rows, index[s + e * (place[j] - place[k])]); push!(cols, col)
                    push!(vals, w)
                end
            end
        end
        push!(rows, col); push!(cols, col); push!(vals, diag)
    end
    A = sparse(rows, cols, vals, length(states), length(states))
    @assert A == transpose(A)
    return A, den
end

"Product states of the ring of `cells` cells grouped by total S^z = M ≥ 0 (M < 0 by spin flip)."
function charge_groups(cells)
    L = 2cells; groups = [Int[] for _ in 0:L]
    for s in 0:3^L-1
        x = s; m = -L
        for _ in 1:L; m += x % 3; x ÷= 3; end
        m >= 0 && push!(groups[m+1], s)
    end
    @assert length(groups[1]) + 2sum(length, groups[2:end]) == 9^cells
    groups
end

"Translation by one cell (two sites) of a base-3 product state."
translate(s, L) = (s ÷ 9) + (s % 9) * 3^(L - 2)
"Site-order reversal of a base-3 product state."
reflect_state(s, L) = foldl((a, _) -> (a[1] * 3 + a[2] % 3, a[2] ÷ 3), 1:L; init = (0, s))[1]

"Real dihedral character basis; generic momenta paired with exact multiplicity two."
function real_orbit_basis(states, cells, k, parity, reflection)
    L = 2cells; flipmax = 3^L - 1; idx = Dict(s => i for (i, s) in enumerate(states))
    visited = Set{Int}(); rows = Int[]; cols = Int[]; vals = Float64[]; col = 0
    generic = !(k == 0 || 2k == cells)
    factor = (generic ? 2 : 1) / (cells * 2 * (parity == 0 ? 1 : 2))
    for s in states
        s in visited && continue
        orbit = Set{Int}()
        for f in 0:(parity == 0 ? 0 : 1), h in 0:1
            t = f == 0 ? s : flipmax - s; t = h == 0 ? t : reflect_state(t, L)
            for r in 0:cells-1; push!(orbit, t); t = translate(t, L); end
        end
        union!(visited, orbit); os = sort!(collect(orbit))
        oi = Dict(t => j for (j, t) in enumerate(os))
        Pr = zeros(length(os), length(os))
        for (j, v) in enumerate(os), f in 0:(parity == 0 ? 0 : 1), h in 0:1
            t = f == 0 ? v : flipmax - v; t = h == 0 ? t : reflect_state(t, L)
            for r in 0:cells-1
                Pr[oi[t], j] += factor * cospi(2 * k * r / cells) * (f == 0 ? 1 : parity) *
                                (h == 0 ? 1 : reflection)
                t = translate(t, L)
            end
        end
        F = eigen(Symmetric(Pr))
        @assert maximum(min.(abs.(F.values), abs.(F.values .- 1))) < 1e-10
        for a in findall(>(0.5), F.values)
            col += 1
            for (i, t) in enumerate(os)
                push!(rows, idx[t]); push!(cols, col); push!(vals, F.vectors[i, a])
            end
        end
    end
    sparse(rows, cols, vals, length(states), col)
end

"Full dense block spectra (floating point); charge and conjugate momenta carry multiplicity."
function dense_spectrum(cells, delta)
    levels = Float64[]; records = Any[]; start = time(); groups = charge_groups(cells)
    for m in 0:2cells
        states = groups[m+1]; A, den = integer_hamiltonian(cells, delta, states)
        for k in 0:fld(cells, 2), p in (m == 0 ? (-1, 1) : (0,)),
            reflection in (k == 0 || 2k == cells ? (-1, 1) : (1,))
            B = real_orbit_basis(states, cells, k, p, reflection); dim = size(B, 2)
            dim == 0 && continue
            H = Matrix(B' * (A * B)) / den
            @assert maximum(abs, H - H') < 1e-10
            mult = (m == 0 ? 1 : 2) * (k == 0 || 2k == cells ? 1 : 2)
            t = @elapsed E = eigvals(Symmetric((H + H') / 2))
            for _ in 1:mult; append!(levels, E); end
            push!(records, (charge = m, momentum = k, spin_flip = p, reflection = reflection,
                dimension = dim, multiplicity = mult, seconds = t))
        end
    end
    @assert sum(r.dimension * r.multiplicity for r in records) == 9^cells == length(levels)
    sort!(levels)
    return PartialSpectrum(levels, 9^cells, true), records, time() - start
end
