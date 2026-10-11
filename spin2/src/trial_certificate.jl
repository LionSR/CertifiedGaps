# Included into Spin2Gap.Trials. All accepting energy comparisons are exact integers.

function quantized_tensors(t, x; scale = 10^8)
    parameters = [round(BigInt, scale*x[i]/t.norms[i]) for i in eachindex(x)]
    return parameters, polynomial_tensors(t, parameters)
end

function polynomial_tensors(t, parameters)
    P = [zeros(BigInt, t.D, t.D) for _ in 1:5]
    for (i, (a, b)) in enumerate(t.pairs), m in 1:5
        P[m][t.ranges[a], t.ranges[b]] .= parameters[i]*t.polynomials[i][m]
    end
    return P
end

function exact_blocks(t, P)
    amax = maximum(maximum(abs, p) for p in P)
    bmax = t.D*amax^2
    @assert 2048bmax^2 < typemax(Int128)
    Q = [Int128.(p) for p in P]
    B = [Q[m]*Q[n] for m in 1:5, n in 1:5]
    charges = last.(t.labels)
    blocks = Any[]
    for delta in (-2maximum(t.spins)):(2maximum(t.spins))
        pairs = [(a, b) for a in 1:t.D for b in 1:t.D if charges[a]-charges[b]==delta]
        T, O =
            zeros(Int128, length(pairs), length(pairs)), zeros(Int128, length(pairs), length(pairs))
        for (row, (a, b)) in enumerate(pairs), (col, (c, d)) in enumerate(pairs)
            T[row, col] = sum(WEIGHTS[m]*Q[m][a, c]*Q[m][b, d] for m in 1:5)
            O[row, col] = sum(
                (m-3)*(n-3)*WEIGHTS[m]*WEIGHTS[n]*B[m, n][a, c]*B[m, n][b, d] for
                m in 1:5, n in 1:5
            )
            for m in -2:1, n in -1:2
                weight = div(WEIGHTS[m+3]*WEIGHTS[n+3]*(2-m)*(2+n), 2)
                O[row, col] +=
                    weight*(B[m+3, n+3][a, c]*B[m+4, n+2][b, d]+B[m+4, n+2][a, c]*B[m+3, n+3][b, d])
            end
        end
        push!(
            blocks,
            (; delta, T = Nemo.matrix(Nemo.ZZ, BigInt.(T)), O = Nemo.matrix(Nemo.ZZ, BigInt.(O))),
        )
    end
    return blocks
end

function exact_contract(blocks, n)
    @assert iseven(n) && n>=4
    U, V = big(0), big(0)
    records = Any[]
    for b in blocks
        X = b.T^(n-2)
        u, v = BigInt(Nemo.tr(X*b.O)), BigInt(Nemo.tr(X*b.T^2))
        U += u
        V += v
        push!(
            records,
            Dict(
                "charge_difference"=>b.delta,
                "dimension"=>Nemo.nrows(b.T),
                "U"=>string(u),
                "V"=>string(v),
            ),
        )
    end
    @assert V>0
    return (; U, V, records)
end

"Direct four-site weighted amplitudes, independently of the double-layer construction."
function physical_integer_check(P)
    words = collect(Iterators.product(ntuple(_->1:5, 4)...))
    amplitudes = Dict(word=>tr(prod(P[m] for m in word)) for word in words)
    V, E = big(0), big(0)
    for (word, value) in amplitudes
        weight = prod(WEIGHTS[m] for m in word)
        V += weight*value^2
        for j in 1:4
            k = mod1(j+1, 4)
            m, n = word[j]-3, word[k]-3
            E += m*n*weight*value^2
            for direction in (-1, 1)
                abs(m+direction)<=2 && abs(n-direction)<=2 || continue
                target = ntuple(i->word[i]+(i==j ? direction : i==k ? -direction : 0), 4)
                factor = direction==1 ? (2-m)*(2+n) : (2+m)*(2-n)
                @assert iseven(weight*factor)
                E += div(weight*factor, 2)*value*amplitudes[target]
            end
        end
    end
    return (; V, E)
end
