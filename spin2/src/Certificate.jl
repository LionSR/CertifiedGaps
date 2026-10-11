module Certificate

using ..Trials, LinearAlgebra, SparseArrays
import Nemo
const ST = Trials

"Exact covariance of the physical two-site energy insertion in polynomial gauge."
function check_physical_scalar()
    plus, minus = ST.ladders(2)
    identity = Matrix{BigInt}(I, 5, 5)
    z = Matrix(Diagonal(big.(-2:2)))
    h2 = 2kron(z, z)+kron(plus, minus)+kron(minus, plus)
    for generator in (plus, minus, z)
        total = kron(generator, identity)+kron(identity, generator)
        @assert total*h2 == h2*total
    end
    metric = Matrix(Diagonal(big.(ST.WEIGHTS)))
    @assert plus*metric == metric*transpose(minus)
    @assert minus*metric == metric*transpose(plus)
    return true
end

"An integer highest-weight matrix in V_j tensor dual(V_k), including its metric."
function highest_matrix(j, k, spin)
    abs(j-k) <= spin <= j+k || throw(ArgumentError("absent representation"))
    jp, jm = ST.ladders(j)
    kp, km = ST.ladders(k)
    pairs = [(r, c) for r in 0:2j for c in 0:2k if r-j-c+k == spin]
    equations = zeros(BigInt, (2j+1)*(2k+1), length(pairs))
    for (column, (r, c)) in enumerate(pairs)
        a = zeros(BigInt, 2j+1, 2k+1)
        a[r+1, c+1] = 1
        equations[:, column] = vec(jp*a-a*kp)
    end
    dim, kernel = Nemo.nullspace(Nemo.matrix(Nemo.QQ, equations))
    @assert dim == 1
    coefficients = [
        ST.frac(BigInt(Nemo.numerator(kernel[i, 1])), BigInt(Nemo.denominator(kernel[i, 1])))
        for i in eachindex(pairs)
    ]
    scale = lcm(denominator.(coefficients))
    R = zeros(BigInt, 2j+1, 2k+1)
    for (i, (r, c)) in enumerate(pairs)
        R[r+1, c+1] = BigInt(scale*coefficients[i])*binomial(big(2k), c)
    end
    R = div.(R, gcd(abs.(vec(R))))
    @assert jp*R-R*transpose(km) == zero(R)
    return R
end

"Check the complete local irreducible basis, not just selected highest vectors."
function representation_basis(j, k)
    jp, jm = ST.ladders(j)
    kp, km = ST.ladders(k)
    columns = Vector{BigInt}[]
    highest = Dict{Int,Matrix{BigInt}}()
    for spin in abs(j-k):(j+k)
        R = highest_matrix(j, k, spin)
        highest[spin] = R
        for weight in spin:-1:(-spin)
            @assert all(iszero(R[r+1, c+1]) || r-j-c+k == weight for r in 0:2j, c in 0:2k)
            push!(columns, vec(R))
            next = jm*R-R*transpose(kp)
            if weight > -spin
                @assert jp*next-next*transpose(km) == (spin+weight)*(spin-weight+1)*R
            else
                @assert iszero(next)
            end
            R = next
        end
    end
    basis = reduce(hcat, columns)
    @assert size(basis) == ((2j+1)*(2k+1), (2j+1)*(2k+1))
    @assert !iszero(Nemo.det(Nemo.matrix(Nemo.ZZ, basis)))
    return highest
end

"Check every tensor covariance entry using the sparse ladder generators."
function check_tensor_covariance(t, P)
    @assert length(P) == 5 && all(size(p) == (t.D, t.D) for p in P)
    previous, following = zeros(Int, t.D), zeros(Int, t.D)
    for range in t.ranges, r in range
        r > first(range) && (previous[r] = r-1)
        r < last(range) && (following[r] = r+1)
    end
    for (i, m) in enumerate(-2:2), r in 1:t.D, c in 1:t.D
        jr, mr = t.labels[r]
        jc, mc = t.labels[c]
        plus =
            (previous[r] == 0 ? big(0) : (jr-mr+1)*P[i][previous[r], c]) -
            (following[c] == 0 ? big(0) : (jc-mc)*P[i][r, following[c]])
        minus =
            (following[r] == 0 ? big(0) : (jr+mr+1)*P[i][following[r], c]) -
            (previous[c] == 0 ? big(0) : (jc+mc)*P[i][r, previous[c]])
        @assert plus == (m == 2 ? 0 : (2-m)*P[i+1][r, c])
        @assert minus == (m == -2 ? 0 : (2+m)*P[i-1][r, c])
    end
    return true
end

function channel_entry_factory(
    Q::Array{T,3},
    B,
    charges,
    ranges,
    diagonal_terms,
    ladder_terms,
) where {T}
    function channel_entry(a, c, b, d, entries)
        tv, ov = zero(T), zero(T)
        @inbounds for (r, s, v) in entries
            rr, ss = ranges[b][r], ranges[d][s]
            delta = charges[a]-charges[rr]
            delta == charges[c]-charges[ss] || continue
            if abs(delta) <= 2
                m = delta+3
                tv += v*ST.WEIGHTS[m]*Q[a, rr, m]*Q[c, ss, m]
            end
            abs(delta) <= 4 || continue
            for (m, n, weight) in diagonal_terms[delta]
                ov += v*weight*B[m, n][a, rr]*B[m, n][c, ss]
            end
            for (m, n, p, q, weight) in ladder_terms[delta]
                ov += v*weight*(B[m, n][a, rr]*B[p, q][c, ss]+B[p, q][a, rr]*B[m, n][c, ss])
            end
        end
        return BigInt(tv), BigInt(ov)
    end
    return channel_entry
end

"Exact multiplicity matrices without constructing the full D^2 transfer matrix."
function reduced_blocks(t, P; verify_images = false, progress = false)
    @assert check_physical_scalar()
    js = unique(t.spins)
    basis = Dict((j, k) => representation_basis(j, k) for j in js, k in js)
    progress && println("complete irreducible bases verified")
    progress && flush(stdout)
    # The tensor intertwining equations are exact. Together with the spin-scalar
    # bond operator they imply covariance of both transfer channels. The checked
    # complete bases above therefore account for every irrep and every weight.
    @assert check_tensor_covariance(t, P)
    amax = maximum(maximum(abs, p) for p in P)
    rsum = maximum(sum(abs, R) for b in values(basis) for R in values(b))
    # First certify native two-tensor products, then use their ACTUAL size to
    # bound the projected sums. The total absolute physical coefficients are
    # 16 for T and 144 + 1024 = 1168 for O.
    product_type = t.D*amax^2 < typemax(Int128) ? Int128 : BigInt
    Q = product_type.(cat(P...; dims = 3))
    sparse_Q = [sparse(Q[:, :, m]) for m in 1:5]
    B = [Matrix(sparse_Q[m]*sparse_Q[n]) for m in 1:5, n in 1:5]
    bmax = BigInt(maximum(maximum(abs, b) for b in B))
    bound = rsum*max(16amax^2, 1168bmax^2)
    integer_type = bound < typemax(Int128) ? Int128 : BigInt
    Q = integer_type.(Q)
    B = [integer_type.(b) for b in B]
    progress && println("exact tensor covariance verified; projected arithmetic=$integer_type")
    progress && flush(stdout)
    charges = last.(t.labels)
    diagonal_terms = Dict(
        delta => [
            (m, n, (m-3)*(n-3)*ST.WEIGHTS[m]*ST.WEIGHTS[n]) for m in 1:5 for
            n in 1:5 if m+n-6 == delta && m != 3 && n != 3
        ] for delta in -4:4
    )
    ladder_terms = Dict(
        delta => [
            (m+3, n+3, m+4, n+2, div(ST.WEIGHTS[m+3]*ST.WEIGHTS[n+3]*(2-m)*(2+n), 2)) for
            m in -2:1 for n in -1:2 if m+n == delta
        ] for delta in -4:4
    )
    channel_entry = channel_entry_factory(Q, B, charges, t.ranges, diagonal_terms, ladder_terms)
    reduced = Any[]
    for spin in 0:2maximum(t.spins)
        pairs = [
            (a, c) for a in eachindex(t.spins) for
            c in eachindex(t.spins) if haskey(basis[(t.spins[a], t.spins[c])], spin)
        ]
        isempty(pairs) && continue
        highest = [basis[(t.spins[a], t.spins[c])][spin] for (a, c) in pairs]
        pivots = [argmax(abs.(R)) for R in highest]
        entries = [
            [
                (r, s, integer_type(R[r, s])) for r in axes(R, 1) for
                s in axes(R, 2) if !iszero(R[r, s])
            ] for R in highest
        ]
        T, O = zeros(ST.F, length(pairs), length(pairs)), zeros(ST.F, length(pairs), length(pairs))
        for (row, (a, c)) in enumerate(pairs), (col, (b, d)) in enumerate(pairs)
            pivot = pivots[row]
            tv, ov = channel_entry(t.ranges[a][pivot[1]], t.ranges[c][pivot[2]], b, d, entries[col])
            T[row, col], O[row, col] = tv//highest[row][pivot], ov//highest[row][pivot]
        end
        if verify_images
            lookup = Dict(pair => i for (i, pair) in enumerate(pairs))
            for (col, (b, d)) in enumerate(pairs),
                a in eachindex(t.spins),
                c in eachindex(t.spins),
                r in 1:t.dims[a],
                s in 1:t.dims[c]

                t.labels[t.ranges[a][r]][2]-t.labels[t.ranges[c][s]][2] == spin || continue
                tv, ov = channel_entry(t.ranges[a][r], t.ranges[c][s], b, d, entries[col])
                row = get(lookup, (a, c), 0)
                if row == 0
                    @assert iszero(tv) && iszero(ov)
                else
                    @assert tv == T[row, col]*highest[row][r, s]
                    @assert ov == O[row, col]*highest[row][r, s]
                end
            end
        end
        push!(
            reduced,
            (; spin, pairs, highest, T = Nemo.matrix(Nemo.QQ, T), O = Nemo.matrix(Nemo.QQ, O)),
        )
        progress && println("exact reduced spin=$spin dimension=$(length(pairs))")
        progress && flush(stdout)
    end
    @assert sum((2b.spin+1)*length(b.pairs) for b in reduced) == t.D^2
    return reduced
end

"Trace of a product without constructing the full matrix product."
function trace_product(A, B)
    @assert Nemo.nrows(A) == Nemo.ncols(A) == Nemo.nrows(B) == Nemo.ncols(B)
    value = zero(Nemo.base_ring(A))
    for i in 1:Nemo.nrows(A), j in 1:Nemo.ncols(A)
        value += A[i, j]*B[j, i]
    end
    return value
end

"Share outward-rounded powers across several power-of-two chain lengths."
function ball_contract_many(reduced, lengths; precision = 192, scale = 1, progress = false)
    @assert !isempty(lengths) && all(n -> n >= 4 && ispow2(n), lengths) && scale > 0
    RR = Nemo.ArbField(precision)
    totals = Dict(n => (U = RR(0), V = RR(0)) for n in lengths)
    for b in reduced
        T = Nemo.matrix(RR, b.T)/RR(big(scale)^2)
        O = Nemo.matrix(RR, b.O)/RR(big(scale)^4)
        T2 = T*T
        X, Y, n = T2, T2*T2, 4
        while n <= maximum(lengths)
            if haskey(totals, n)
                old = totals[n]
                totals[n] = (
                    U = old.U+(2b.spin+1)*trace_product(X, O),
                    V = old.V+(2b.spin+1)*trace_product(X, T2),
                )
            end
            n == maximum(lengths) && break
            X, Y = X*Y, Y*Y
            n *= 2
        end
        progress && println("ball contraction spin=$(b.spin) dimension=$(length(b.pairs))")
        progress && flush(stdout)
    end
    @assert all(c.V > 0 for c in values(totals))
    return Dict(n => (; c.U, c.V, energy = c.U/c.V) for (n, c) in totals)
end

"Split the exact commuting involution X -> U X^T U^T, U a pi rotation."
function split_transpose_rotation(reduced; progress = false)
    result = Any[]
    for b in reduced
        lookup = Dict(pair=>i for (i, pair) in enumerate(b.pairs))
        partner = [lookup[(c, a)] for (a, c) in b.pairs]
        signs = Int[]
        for (i, R) in enumerate(b.highest)
            transformed = [(-1)^(r+c)*R[end-c+1, end-r+1] for r in axes(R, 2), c in axes(R, 1)]
            target = b.highest[partner[i]]
            pivot = argmax(abs.(target))
            sign = transformed[pivot]//target[pivot]
            @assert sign in (-1, 1) && transformed == sign*target
            push!(signs, Int(sign))
        end
        @assert all(signs[i]*signs[partner[i]] == 1 for i in eachindex(partner))
        # Direct exact commutation checks make the matrix split independent of
        # floating-point gauges and of an assumed phase convention for the bases.
        for M in (b.T, b.O), i in eachindex(partner), j in eachindex(partner)
            @assert signs[partner[i]]*M[partner[i], j] == signs[j]*M[i, partner[j]]
        end
        for parity in (1, -1)
            representatives = [
                i for i in eachindex(partner) if
                i < partner[i] || (i == partner[i] && signs[i] == parity)
            ]
            isempty(representatives) && continue
            function restrict(M)
                return Nemo.matrix(
                    Nemo.QQ,
                    [
                        M[i, j]+(j == partner[j] ? 0 : parity*signs[j]*M[i, partner[j]]) for
                        i in representatives, j in representatives
                    ],
                )
            end
            push!(
                result,
                (;
                    spin = b.spin,
                    parity,
                    pairs = b.pairs[representatives],
                    T = restrict(b.T),
                    O = restrict(b.O),
                ),
            )
        end
        progress && println("transpose-rotation split spin=$(b.spin) dimension=$(length(b.pairs))")
        progress && flush(stdout)
    end
    @assert sum((2b.spin+1)*length(b.pairs) for b in result) ==
            sum((2b.spin+1)*length(b.pairs) for b in reduced)
    return result
end

function exact_contract(reduced, n)
    @assert iseven(n) && n >= 4
    U, V = zero(Nemo.QQ), zero(Nemo.QQ)
    for b in reduced
        X = b.T^(n-2)
        U += (2b.spin+1)*Nemo.tr(X*b.O)
        V += (2b.spin+1)*Nemo.tr(X*b.T^2)
    end
    @assert Nemo.denominator(U) == 1 && Nemo.denominator(V) == 1
    return (; U = BigInt(Nemo.numerator(U)), V = BigInt(Nemo.numerator(V)))
end

"Outward-rounded finite variational enclosure for an exact integer witness."
function ball_contract(reduced, n; precision = 192, scale = 1)
    @assert iseven(n) && n >= 4 && scale > 0
    RR = Nemo.ArbField(precision)
    U, V = RR(0), RR(0)
    for b in reduced
        T = Nemo.matrix(RR, b.T)/RR(big(scale)^2)
        O = Nemo.matrix(RR, b.O)/RR(big(scale)^4)
        X = T^(n-2)
        U += (2b.spin+1)*Nemo.tr(X*O)
        V += (2b.spin+1)*Nemo.tr(X*T^2)
    end
    @assert V > 0
    return (; U, V, energy = U/V)
end

end
