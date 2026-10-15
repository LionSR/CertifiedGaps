# Complete symmetry-resolved spectra of the ladder ring in the original spin basis.
#
# Bits 2j and 2j+1 of a basis state are the two spins of rung j (bit value one is spin
# down). The basis is resolved by down-spin number m, rung momentum k, leg exchange p and,
# at m = L, global spin inversion q. Orbits and stabilisers are exact permutations; a
# character is kept exactly when it is trivial on the stabiliser.

"Translate a basis state by one rung on a ring of `L` rungs."
translate(s, L) = ((s << 2) & ((1 << (2L)) - 1)) | (s >> (2L - 2))

"Exchange the two legs of every rung."
function swaplegs(s, L)
    mask = sum(1 << (2j) for j in 0:L-1)
    ((s & mask) << 1) | ((s >> 1) & mask)
end

"Reflect the ring of rungs, j -> -j mod L."
function reflectrungs(s, L)
    out = 0
    for j in 0:L-1
        out |= ((s >> (2j)) & 3) << (2mod(-j, L))
    end
    out
end

"""
    orbit_data(L, m)

Exact orbits of the states with `m` down spins under C_L × leg exchange (× spin reversal at
m = L): representatives, orbit sizes, stabilisers as `(t, a, b)` group elements, and for
every state its representative index `rid` and the group element `(tt, aa, bb)` reaching it.
"""
function orbit_data(L, m)
    N = 1 << (2L)
    rid = zeros(Int, N); tt = zeros(Int, N); aa = zeros(Int, N); bb = zeros(Int, N)
    reps = Int[]; sizes = Int[]; stabs = Vector{NTuple{3,Int}}[]
    flip = m == L
    for s in 0:N-1
        count_ones(s) == m && rid[s+1] == 0 || continue
        push!(reps, s); r = length(reps); o = 0; stab = NTuple{3,Int}[]
        for b in 0:Int(flip), a in 0:1
            v = b == 0 ? s : xor(s, N - 1)
            v = a == 0 ? v : swaplegs(v, L)
            for t in 0:L-1
                v == s && push!(stab, (t, a, b))
                if rid[v+1] == 0
                    rid[v+1] = r; tt[v+1] = t; aa[v+1] = a; bb[v+1] = b; o += 1
                end
                v = translate(v, L)
            end
        end
        push!(sizes, o); push!(stabs, stab)
    end
    @assert sum(sizes) == binomial(2L, m)
    (; L, m, reps, sizes, stabs, rid, tt, aa, bb, flip)
end

"Character exponent (in units of π/L) of group element (t, a, b) in sector (k, p, q)."
phase_exponent(k, p, q, t, a, b, L) = mod(2k * t + L * (p * a + q * b), 2L)

"Representatives whose stabiliser is trivial under the (k, p, q) character, and their index."
function block_basis(o, k, p, q)
    selected = findall(st -> all(phase_exponent(k, p, q, t, a, b, o.L) == 0
                                 for (t, a, b) in st), o.stabs)
    index = zeros(Int, length(o.reps))
    for (i, r) in enumerate(selected); index[r] = i; end
    selected, index
end

"Exact Heisenberg matrix elements ⟨t|H|s⟩ in the original spin basis, as rationals."
function transitions(s, L, J)
    out = Dict{Int,Rational{Int}}(s => 0 // 1)
    for j in 0:L-1
        pairs = ((2j, 2j + 1, J), (2j, 2mod(j + 1, L), 1 // 1),
                 (2j + 1, 2mod(j + 1, L) + 1, 1 // 1))
        for (a, b, w) in pairs
            different = ((s >> a) & 1) != ((s >> b) & 1)
            out[s] += different ? -w / 4 : w / 4
            if different
                t = xor(s, (1 << a) | (1 << b)); out[t] = get(out, t, 0 // 1) + w / 2
            end
        end
    end
    out
end

"Multiplicity of the (m, k) block: opposite magnetisations and conjugate momenta coincide."
block_multiplicity(m, k, L) = (m == L ? 1 : 2) * (mod(2k, L) == 0 ? 1 : 2)

"""
    verified_block(o, k, p, q, J) -> (H, herror)

Float64 centre `H` of the exact Hermitian symmetry block, with `herror` an upper bound on
‖H_exact - H‖₂. Entries are rational couplings times enclosed roots of unity and square
roots of exact orbit-size ratios, accumulated in interval arithmetic; `herror` is the
largest row sum of the entrywise radii, which bounds the 2-norm of a Hermitian error.
"""
function verified_block(o, k, p, q, J)
    sel, idx = block_basis(o, k, p, q); d = length(sel)
    realblock = mod(2k, o.L) == 0
    H = zeros(realblock ? Float64 : ComplexF64, d, d)
    re = fill(interval(0.0), d, d); im = copy(re)
    phases = [(cospi(to_interval(e // o.L)), -sinpi(to_interval(e // o.L))) for e in 0:2o.L-1]
    for (col, r) in enumerate(sel)
        for (s, h) in transitions(o.reps[r], o.L, J)
            rr = o.rid[s+1]; row = idx[rr]; row == 0 && continue
            e = phase_exponent(k, p, q, o.tt[s+1], o.aa[s+1], o.bb[s+1], o.L)
            amp = to_interval(h) * sqrt(to_interval(o.sizes[r] // o.sizes[rr]))
            re[row, col] += amp * phases[e+1][1]
            im[row, col] += amp * phases[e+1][2]
        end
    end
    errors = fill(interval(0.0), d)
    for j in 1:d, i in 1:j
        a = mid(re[i, j]); b = realblock || i == j ? 0.0 : mid(im[i, j])
        H[i, j] = realblock ? a : complex(a, b); H[j, i] = conj(H[i, j])
        err = abs(re[i, j] - a) + abs(im[i, j] - b)
        # The exact block is Hermitian, so the upper triangle determines both halves.
        errors[i] += err
        i != j && (errors[j] += err)
    end
    H, sup(interval(0.0, maximum(sup, errors; init = 0.0)))
end

"""
    numerical_block(o, k, p, q, J) -> Matrix{Float64}

Unverified real symmetric form of a symmetry block, for floating-point scans only. Complex
momentum blocks are made real by the antiunitary rung reflection times complex conjugation;
the change of basis is sparse to avoid large dense complex intermediates at ten rungs.
"""
function numerical_block(o, k, p, q, J)
    sel, idx = block_basis(o, k, p, q); d = length(sel)
    ii = Int[]; jj = Int[]; vv = ComplexF64[]
    phases = [cispi(-e / o.L) for e in 0:2o.L-1]
    for (col, r) in enumerate(sel), (s, h) in transitions(o.reps[r], o.L, J)
        rr = o.rid[s+1]; row = idx[rr]; row == 0 && continue
        e = phase_exponent(k, p, q, o.tt[s+1], o.aa[s+1], o.bb[s+1], o.L)
        push!(ii, row); push!(jj, col)
        push!(vv, Float64(h) * sqrt(o.sizes[r] / o.sizes[rr]) * phases[e+1])
    end
    H = sparse(ii, jj, vv, d, d)
    @assert isapprox(H, H'; atol = 1e-10, rtol = 1e-12)
    H = sparse(Hermitian(H))
    mod(2k, o.L) == 0 && return Matrix(real(H))
    ii = Int[]; jj = Int[]; vv = ComplexF64[]; seen = falses(d)
    for (i, r) in enumerate(sel)
        seen[i] && continue
        s = reflectrungs(o.reps[r], o.L); j = idx[o.rid[s+1]]
        @assert j != 0
        e = phase_exponent(k, p, q, o.tt[s+1], o.aa[s+1], o.bb[s+1], o.L)
        t = cispi(-e / (2o.L)); seen[i] = true; seen[j] = true
        if i == j
            push!(ii, i); push!(jj, i); push!(vv, t)
        else
            append!(ii, (i, j, i, j)); append!(jj, (i, i, j, j))
            append!(vv, (t / sqrt(2), t / sqrt(2), im * t / sqrt(2), -im * t / sqrt(2)))
        end
    end
    U = sparse(ii, jj, vv, d, d); A = U' * H * U
    @assert maximum(abs, imag(A); init = 0.0) < 1e-9
    Matrix(real(A))
end

"Upper bound on the Frobenius norm, accumulated in interval arithmetic from exact entries."
function frob(A)
    s = interval(0.0)
    for z in A
        s += interval(real(z))^2 + interval(imag(z))^2
    end
    sqrt(s)
end

"""
    verified_eigen(H, herror) -> (λ, ρ, audit)

Eigenvalues `λ` of a LAPACK eigendecomposition of `H` and a radius `ρ` such that the sorted
eigenvalues of every Hermitian matrix within `herror` of `H` lie within `ρ` of `λ`. With
residual r ≥ ‖HQ - QΛ‖ and η ≥ ‖Q*Q - I‖ < 1, the polar factor U of Q obeys
‖Q - U‖ ≤ η / (1 + √(1-η)), and Weyl gives ρ = r + (‖H‖ + max|λ|)‖Q - U‖ + herror.
"""
function verified_eigen(H, herror)
    n = size(H, 1); f = eigen(Hermitian(H)); V = f.vectors; lam = f.values
    R = H * V - V * Diagonal(lam); E = V' * V - I
    # R and E are computed in floating point: γ_{8n+32} bounds the real operations in
    # complex dot products plus the subtraction, for any summation order.
    u = interval(2.0)^(-53); g = (8n + 32) * u / (1 - (8n + 32) * u)
    hn = frob(H); qn = frob(V); ln = frob(lam)
    # Absolute underflow allowance, far larger than 8n^3 smallest subnormals.
    tiny = interval(2.0)^(-1000)
    r = frob(R) + g * (hn * qn + qn * ln) + tiny
    eta = frob(E) + g * (qn^2 + sqrt(interval(n))) + tiny
    @assert sup(eta) < 1
    dist = eta / (1 + sqrt(1 - eta))
    rho = r + (hn + maximum(abs, lam)) * dist + interval(herror)
    audit = Dict("dimension" => n, "matrix_error" => string(Rational{BigInt}(herror)),
        "residual_upper" => string(Rational{BigInt}(sup(r))),
        "orthogonality_upper" => string(Rational{BigInt}(sup(eta))),
        "eigenvalue_radius" => string(Rational{BigInt}(sup(rho))))
    lam, sup(rho), audit
end

"""
    full_spectrum(L, J; verified = false)

All 4^L eigenvalues of the ring of `L` rungs, sorted, with exact block multiplicities. With
`verified = true` every eigenvalue carries a rigorous radius (`radii`); otherwise radii are
zero and the result is floating-point evidence only. Returns
`(; energies, radii, rows, seconds, dimension)`, `rows` holding per-block audit records.
"""
function full_spectrum(L, J; verified = false)
    start = time(); energies = Float64[]; radii = Float64[]; rows = Any[]; dimension = 0
    for m in 0:L
        o = orbit_data(L, m)
        for k in 0:fld(L, 2), p in 0:1, q in 0:Int(o.flip)
            mult = block_multiplicity(m, k, L)
            tick = time()
            if verified
                H, err = verified_block(o, k, p, q, J)
                d = size(H, 1); d == 0 && continue
                ev, rho, audit = verified_eigen(H, err)
            else
                H = numerical_block(o, k, p, q, J)
                d = size(H, 1); d == 0 && continue
                ev = eigvals(Symmetric(H)); rho = 0.0; audit = Dict{String,Any}()
            end
            merge!(audit, Dict("m" => m, "k" => k, "leg_parity" => (-1)^p,
                "spin_parity" => o.flip ? (-1)^q : 0, "dimension" => d,
                "multiplicity" => mult, "seconds" => time() - tick))
            push!(rows, audit)
            for _ in 1:mult; append!(energies, ev); append!(radii, fill(rho, d)); end
            dimension += mult * d
        end
        @printf("J=%g L=%d m=%d levels=%d elapsed=%.1fs\n", Float64(J), L, m, dimension,
            time() - start)
        flush(stdout)
    end
    @assert dimension == 4^L == length(energies)
    perm = sortperm(energies); energies = energies[perm]; radii = radii[perm]
    (; energies, radii, rows, seconds = time() - start, dimension)
end

"""
    sector_dimensions(L) -> Dict((m, k, leg_parity, spin_parity) => (dimension, multiplicity))

Exact block dimensions and multiplicities, independent of any eigensolve or saved file.
"""
function sector_dimensions(L)
    rows = Dict{NTuple{4,Int},Tuple{Int,Int}}()
    for m in 0:L
        o = orbit_data(L, m)
        for k in 0:fld(L, 2), p in 0:1, r in 0:Int(o.flip)
            selected, _ = block_basis(o, k, p, r)
            isempty(selected) && continue
            rows[(m, k, (-1)^p, o.flip ? (-1)^r : 0)] =
                (length(selected), block_multiplicity(m, k, L))
        end
    end
    @assert sum(d * m for (d, m) in values(rows)) == 4^L
    return rows
end

"""
    enclosed_logZ(spectrum, β) -> Interval

Enclosure of log Z_L(β) = log Σ exp(-β E) over every eigenvalue enclosure `E ± r` of a
verified spectrum, shifted by the lowest centre to avoid overflow.
"""
function enclosed_logZ(s, beta)
    shift = s.energies[1]
    z = interval(0.0)
    for (e, r) in zip(s.energies, s.radii)
        z += exp(-to_interval(beta) * ((interval(e) + interval(-r, r)) - shift))
    end
    log(z) - to_interval(beta) * shift
end
