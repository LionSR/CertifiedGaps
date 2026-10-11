module Trials

using LinearAlgebra, Random, SHA, Dates, JSON3
import Nemo
using ThermalMPOGap.Spins: integer_ladders, binomial_weights
const F = Rational{BigInt}
frac(a, b = 1) = BigInt(a)//BigInt(b)
const WEIGHTS = binomial_weights(2) # binomial(4,2+m), m=-2:2, i.e. [1, 4, 6, 4, 1]

"Integer ladders in the binomial gauge, m = -j:j ascending (`ThermalMPOGap.Spins`)."
const ladders = integer_ladders

"Integer rank-two spherical tensor from exact highest-weight equations."
function intertwiner(j, k)
    abs(j-k) <= 2 <= j+k || throw(ArgumentError("forbidden spin-two coupling"))
    jp, jm = ladders(j)
    kp, km = ladders(k)
    pairs = [(r, c) for r in 0:2j for c in 0:2k if (r-j)-(c-k)==2]
    equations = zeros(BigInt, (2j+1)*(2k+1), length(pairs))
    for (col, (r, c)) in enumerate(pairs)
        M = zeros(BigInt, 2j+1, 2k+1)
        M[r+1, c+1] = 1
        equations[:, col] = vec(jp*M-M*kp)
    end
    dim, N = Nemo.nullspace(Nemo.matrix(Nemo.QQ, equations))
    @assert dim == 1
    highest = zeros(F, 2j+1, 2k+1)
    for (i, (r, c)) in enumerate(pairs)
        highest[r+1, c+1] = frac(BigInt(Nemo.numerator(N[i, 1])), BigInt(Nemo.denominator(N[i, 1])))
    end
    P = [zeros(F, 2j+1, 2k+1) for _ in -2:2]
    P[5] = highest
    for m in 2:-1:-1
        P[m+2] = (jm*P[m+3]-P[m+3]*km)/(2+m)
    end
    scale = lcm(denominator.(vcat(vec.(P)...)))
    ints = [BigInt.(p*scale) for p in P]
    divisor = gcd(abs.(vcat(vec.(ints)...)))
    ints = [div.(p, divisor) for p in ints]
    for (i, m) in enumerate(-2:2)
        @assert jp*ints[i]-ints[i]*kp == (m==2 ? zero(ints[i]) : (2-m)*ints[i+1])
        @assert jm*ints[i]-ints[i]*km == (m == -2 ? zero(ints[i]) : (2+m)*ints[i-1])
    end
    return ints
end

"Integer virtual multiplets, with multiplicities; exact templates and numerical controls."
function templates(spins::Vector{Int})
    @assert all(>=(0), spins)
    dims = 2spins .+ 1
    ends = cumsum(dims)
    ranges = [(ends[i]-dims[i]+1):ends[i] for i in eachindex(spins)]
    labels = [(j, m) for j in spins for m in (-j):j]
    D = sum(dims)
    pairs, polynomials, numerical, norms = Tuple{Int,Int}[], Any[], Any[], Float64[]
    for (a, j) in enumerate(spins), (b, k) in enumerate(spins)
        abs(j-k) <= 2 <= j+k || continue
        P = intertwiner(j, k)
        # Similarity to an orthonormal virtual basis; physical factors sqrt(w_m).
        left = [1/sqrt(binomial(2j, r)) for r in 0:2j]
        right = [1/sqrt(binomial(2k, r)) for r in 0:2k]
        A = [sqrt(WEIGHTS[m])*Float64.(P[m]) .* (left ./ transpose(right)) for m in 1:5]
        norm = sqrt(sum(sum(abs2, a) for a in A))
        push!(pairs, (a, b))
        push!(polynomials, P)
        push!(numerical, [a/norm for a in A])
        push!(norms, norm)
    end
    singlets = [(a, c) for a in eachindex(spins), c in eachindex(spins) if spins[a]==spins[c]]
    lookup = Dict(pair=>i for (i, pair) in enumerate(pairs))
    entries = [
        (row, col, lookup[(a, b)], lookup[(c, d)], 1/dims[a]) for
        (row, (a, c)) in enumerate(singlets), (col, (b, d)) in enumerate(singlets) if
        haskey(lookup, (a, b)) && haskey(lookup, (c, d))
    ]
    return (;
        spins,
        dims,
        ranges,
        labels,
        D,
        pairs,
        polynomials,
        numerical,
        norms,
        singlets,
        entries,
    )
end

function tensors(t, x)
    @assert length(x)==length(t.pairs)
    A = [zeros(t.D, t.D) for _ in 1:5]
    for (i, (a, b)) in enumerate(t.pairs), m in 1:5
        A[m][t.ranges[a], t.ranges[b]] .= x[i]*t.numerical[i][m]
    end
    return A
end

"Numerical infinite-MPS variational energy; not a gap or an exact certificate."
function singlet_transfer(t, x)
    # SU(2)-invariant fixed points reduce to scalar identities in each multiplet.
    W = zeros(length(t.singlets), length(t.singlets))
    for (row, col, a, b, factor) in t.entries
        W[row, col] = factor*x[a]*x[b]
    end
    return W
end

function canonicalize(t, x)
    W = singlet_transfer(t, x)
    ev = eigen(W)
    i = argmax(real.(ev.values))
    r = real.(ev.vectors[:, i])
    R = zeros(length(t.spins), length(t.spins))
    for (k, (a, b)) in enumerate(t.singlets)
        R[a, b] = r[k]
    end
    R = (R+transpose(R))/2
    tr(R)<0 && (R = -R)
    root, inverse = zero(R), zero(R)
    for j in unique(t.spins)
        inds = findall(==(j), t.spins)
        s = eigen(Symmetric(R[inds, inds]))
        minimum(s.values)>1e-14 || return x
        root[inds, inds] = s.vectors*Diagonal(sqrt.(s.values))*transpose(s.vectors)
        inverse[inds, inds] = s.vectors*Diagonal(inv.(sqrt.(s.values)))*transpose(s.vectors)
    end
    X = zeros(length(t.spins), length(t.spins))
    for (i, (a, b)) in enumerate(t.pairs)
        X[a, b] = x[i]
    end
    Y = inverse*X*root
    result = [Y[a, b] for (a, b) in t.pairs]
    result /= norm(result)
    @assert abs(energy(t, result)-energy(t, x))<1e-9
    return result
end

function energy(t, x)
    W = singlet_transfer(t, x)
    er, el = eigen(W), eigen(transpose(W))
    ir, il = argmax(real.(er.values)), argmax(real.(el.values))
    lambda = real(er.values[ir])
    lambda > 0 || return Inf
    r, l = real.(er.vectors[:, ir]), real.(el.vectors[:, il])
    dot(l, r) != 0 || return Inf
    l /= dot(l, r)
    R, L = zeros(t.D, t.D), zeros(t.D, t.D)
    for (i, (a, b)) in enumerate(t.singlets), k in 1:t.dims[a]
        R[t.ranges[a][k], t.ranges[b][k]] = r[i]
        L[t.ranges[a][k], t.ranges[b][k]] = l[i]/t.dims[a]
    end
    A = tensors(t, x)
    B = [A[m]*A[n] for m in 1:5, n in 1:5]
    LB = [L*B[m, n]*transpose(R) for m in 1:5, n in 1:5]
    value = sum((m-3)*(n-3)*dot(LB[m, n], B[m, n]) for m in 1:5, n in 1:5)
    for m in -2:1, n in -1:2
        value += sqrt((6-m*(m+1))*(6-n*(n-1)))*dot(LB[m+3, n+3], B[m+4, n+2])
    end
    return value/lambda^2
end

"Finite-difference BFGS refinement of numerical candidates only."
function refine(t, x; iterations = 250)
    x = x/norm(x)
    f = y -> energy(t, y)
    function gradient(y)
        g = similar(y)
        for i in eachindex(y)
            step = 1e-5
            yp, ym = copy(y), copy(y)
            yp[i] += step
            ym[i] -= step
            g[i] = (f(yp)-f(ym))/(2step)
        end
        return g
    end
    d = length(x)
    H = Matrix{Float64}(I, d, d)
    value, g = f(x), gradient(x)
    completed = 0
    for iteration in 1:iterations
        norm(g, Inf) < 2e-7 && break
        direction = -H*g
        dot(g, direction) < 0 || (direction = -g; H = Matrix{Float64}(I, d, d))
        norm(direction)>1 && (direction /= norm(direction))
        step = 1.0
        candidate = x+step*direction
        fv = f(candidate)
        while fv > value+1e-4*step*dot(g, direction) && step>1e-9
            step /= 2
            candidate = x+step*direction
            fv = f(candidate)
        end
        step <= 1e-9 && break
        gn = gradient(candidate)
        s, y = candidate-x, gn-g
        sy = dot(s, y)
        if sy > 1e-14
            rho = inv(sy)
            Q = Matrix{Float64}(I, d, d)-rho*s*transpose(y)
            H = Q*H*transpose(Q)+rho*s*transpose(s)
        end
        x, value, g = candidate, fv, gn
        completed = iteration
    end
    x /= norm(x)
    return (; x, energy = f(x), iterations = completed, gradient_norm = norm(g, Inf))
end

"Deterministic Nelder-Mead search; gauge/normalization freedom is left explicit."
function search(t; seed = 1, iterations = 2500, start = nothing)
    rng = MersenneTwister(seed)
    d = length(t.pairs)
    x0 = isnothing(start) ? randn(rng, d) : copy(start)
    x0 /= norm(x0)
    eye = Matrix{Float64}(I, d, d)
    simplex = hcat(x0, [x0+0.12*eye[:, i] for i in 1:d]...)
    values = [energy(t, x) for x in eachcol(simplex)]
    for iteration in 1:iterations
        order = sortperm(values)
        simplex, values = simplex[:, order], values[order]
        center = vec(sum(simplex[:, 1:d]; dims = 2))/d
        reflection = 2center-simplex[:, end]
        fr = energy(t, reflection)
        if fr < values[1]
            expansion = center+2*(reflection-center)
            fe = energy(t, expansion)
            simplex[:, end], values[end] = fe<fr ? (expansion, fe) : (reflection, fr)
        elseif fr < values[d]
            simplex[:, end], values[end] = reflection, fr
        else
            outside = fr < values[end]
            contraction = center+0.5*((outside ? reflection : simplex[:, end])-center)
            fc = energy(t, contraction)
            if fc < (outside ? fr : values[end])
                simplex[:, end], values[end] = contraction, fc
            else
                for i in 2:(d+1)
                    simplex[:, i] = (simplex[:, i]+simplex[:, 1])/2
                    values[i] = energy(t, simplex[:, i])
                end
            end
        end
        maximum(values)-minimum(values) < 1e-12 && iteration>100 && break
    end
    best = argmin(values)
    x = simplex[:, best]/norm(simplex[:, best])
    return (; x, energy = energy(t, x))
end

include("trial_certificate.jl")

end
