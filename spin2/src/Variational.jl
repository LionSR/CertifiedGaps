module Variational

using ..Trials, LinearAlgebra, Random
const ST = Trials

"Numerical angular coefficient for the two-site energy channel."
function angular_coefficient(j, e, k, f, tensor_cache)
    function local_tensors(a, b)
        get!(tensor_cache, (a, b)) do
            P = ST.intertwiner(a, b)
            left = [1 / sqrt(binomial(2a, r)) for r in 0:2a]
            right = [1 / sqrt(binomial(2b, r)) for r in 0:2b]
            A = [sqrt(ST.WEIGHTS[m]) * Float64.(P[m]) .* (left ./ transpose(right)) for m in 1:5]
            factor = sqrt(sum(sum(abs2, a) for a in A))
            [a / factor for a in A]
        end
    end
    ae, eb = local_tensors(j, e), local_tensors(e, k)
    af, fb = local_tensors(j, f), local_tensors(f, k)
    B = [ae[m] * eb[n] for m in 1:5, n in 1:5]
    C = [af[m] * fb[n] for m in 1:5, n in 1:5]
    value = sum((m - 3) * (n - 3) * dot(B[m, n], C[m, n]) for m in 1:5, n in 1:5)
    for m in -2:1, n in -1:2
        factor = sqrt((6 - m * (m + 1)) * (6 - n * (n - 1))) / 2
        value += factor * (dot(B[m+3, n+3], C[m+4, n+2]) + dot(B[m+4, n+2], C[m+3, n+3]))
    end
    return value / (2j + 1)
end

"Precompute the quartic energy channel on SU(2) multiplicity spaces."
function polynomial_model(spins::Vector{Int})
    t = ST.templates(spins)
    lookup = Dict(pair => i for (i, pair) in enumerate(t.pairs))
    paths = Dict(
        (a, b) => [
            (e, lookup[(a, e)], lookup[(e, b)]) for
            e in eachindex(spins) if haskey(lookup, (a, e)) && haskey(lookup, (e, b))
        ] for a in eachindex(spins), b in eachindex(spins)
    )
    angular = Dict{NTuple{4,Int},Float64}()
    tensor_cache = Dict{Tuple{Int,Int},Vector{Matrix{Float64}}}()
    quartic = Tuple{Int,Int,Int,Int,Int,Int,Float64}[]
    for (row, (a, c)) in enumerate(t.singlets),
        (col, (b, d)) in enumerate(t.singlets),
        (e, i1, i2) in paths[(a, b)],
        (f, i3, i4) in paths[(c, d)]

        key = (spins[a], spins[e], spins[b], spins[f])
        coefficient = get!(angular, key) do
            angular_coefficient(key..., tensor_cache)
        end
        coefficient == 0 && continue
        push!(quartic, (row, col, i1, i2, i3, i4, coefficient))
    end
    return (; t, quartic, angular)
end

"Factor the quartic channel into products of small multiplicity matrices."
function model(spins::Vector{Int})
    t = ST.templates(spins)
    js = sort(unique(spins))
    groups = Dict(j => findall(==(j), spins) for j in js)
    lookup = Dict(pair => i for (i, pair) in enumerate(t.pairs))
    indices = Dict(
        (j, k) => [lookup[(a, b)] for a in groups[j], b in groups[k]] for
        j in js, k in js if abs(j-k) <= 2 <= j+k
    )
    singlet_lookup = Dict(pair => i for (i, pair) in enumerate(t.singlets))
    blocks = Dict(j => [singlet_lookup[(a, c)] for a in groups[j] for c in groups[j]] for j in js)
    tensor_cache = Dict{Tuple{Int,Int},Vector{Matrix{Float64}}}()
    couplings = Any[]
    equivalent_terms = 0
    for j in js, k in js
        middle = [e for e in js if haskey(indices, (j, e)) && haskey(indices, (e, k))]
        isempty(middle) && continue
        C = [angular_coefficient(j, e, k, f, tensor_cache) for e in middle, f in middle]
        @assert norm(C-transpose(C), Inf) < 1e-12
        push!(couplings, (; j, k, middle, C))
        equivalent_terms +=
            length(groups[j])^2 * length(groups[k])^2 * sum(length(groups[e]) for e in middle)^2
    end
    symmetric_columns = [
        (i, singlet_lookup[(c, a)], a == c ? 1.0 : inv(sqrt(2.0))) for
        (i, (a, c)) in enumerate(t.singlets) if a <= c
    ]
    return (; t, js, groups, indices, blocks, couplings, equivalent_terms, symmetric_columns)
end

"Restrict a transpose-preserving real channel to symmetric multiplicity matrices."
function symmetric_restriction(A, columns)
    result = zeros(length(columns), length(columns))
    @inbounds for (i, (a, b, u)) in enumerate(columns), (j, (c, d, v)) in enumerate(columns)
        result[i, j] =
            u*v*(
                A[a, c]+(a == b ? 0 : A[b, c])+(c == d ? 0 : A[a, d])+(
                    a == b || c == d ? 0 : A[b, d]
                )
            )
    end
    return result
end

function symmetric_expansion(v, columns, dimension)
    result = zeros(dimension)
    @inbounds for (i, (a, b, u)) in enumerate(columns)
        result[a] = u*v[i]
        a == b || (result[b] = u*v[i])
    end
    return result
end

function fixed_points(W)
    er = eigen(W)
    ir = argmax(real.(er.values))
    lambda = real(er.values[ir])
    lambda > 0 || throw(ArgumentError("nonpositive transfer radius"))
    r = real.(er.vectors[:, ir])
    # r'*(lambda*I-W')=0, so this rank-one stabilized solve enforces r'l=1.
    l = (lambda*I-transpose(W)+r*transpose(r))\r
    abs(dot(l, r)) > 1e-12 || throw(ArgumentError("ill-conditioned transfer fixed point"))
    l /= dot(l, r)
    return lambda, r, l
end

function factorized_channels(model, x)
    X = Dict(key => x[indices] for (key, indices) in model.indices)
    W = ST.singlet_transfer(model.t, x)
    B = zero(W)
    paths = Any[]
    for coupling in model.couplings
        (; j, k, middle, C) = coupling
        Y = [X[(j, e)] * X[(e, k)] for e in middle]
        block = zeros(length(model.blocks[j]), length(model.blocks[k]))
        for e in eachindex(middle), f in eachindex(middle)
            block .+= C[e, f] .* kron(Y[e], Y[f])
        end
        B[model.blocks[j], model.blocks[k]] = block
        push!(paths, Y)
    end
    return W, B, X, paths
end

function channels(model, x)
    if hasproperty(model, :couplings)
        W, B, _, _ = factorized_channels(model, x)
        return W, B
    end
    W = ST.singlet_transfer(model.t, x)
    B = zero(W)
    @inbounds for (row, col, i, j, k, l, c) in model.quartic
        B[row, col] += c * x[i] * x[j] * x[k] * x[l]
    end
    return W, B
end

"Analytic derivative including both transfer fixed points; numerical, not certified."
function energy_gradient(model, x; symmetric = true)
    factorized = hasproperty(model, :couplings)
    if factorized
        W, B, X, paths = factorized_channels(model, x)
    else
        W, B = channels(model, x)
    end
    reduced = symmetric && hasproperty(model, :symmetric_columns)
    if reduced
        W = symmetric_restriction(W, model.symmetric_columns)
        B = symmetric_restriction(B, model.symmetric_columns)
    end
    lambda, r, l = fixed_points(W)
    projector = r * transpose(l)
    resolvent_factor = lu(lambda*I-W+projector)
    value = dot(l, B * r) / lambda^2
    br, bl = B*r, transpose(B)*l
    right_response = resolvent_factor\br-projector*br
    left_response = transpose(resolvent_factor)\bl-transpose(projector)*bl
    if reduced
        expand(v) = symmetric_expansion(v, model.symmetric_columns, length(model.t.singlets))
        r, l, right_response, left_response =
            expand(r), expand(l), expand(right_response), expand(left_response)
    end
    gradient = zeros(length(x))
    if factorized
        gx = Dict(key => zero(value) for (key, value) in X)
        for (coupling, Y) in zip(model.couplings, paths)
            (; j, k, middle, C) = coupling
            # singlet vectors enumerate the first matrix index outside the second.
            L = transpose(
                reshape(l[model.blocks[j]], length(model.groups[j]), length(model.groups[j])),
            )
            R = transpose(
                reshape(r[model.blocks[k]], length(model.groups[k]), length(model.groups[k])),
            )
            for e in eachindex(middle)
                gy = zero(Y[e])
                for f in eachindex(middle)
                    gy .+=
                        (C[e, f] * L * Y[f] * transpose(R) + C[f, e] * transpose(L) * Y[f] * R) /
                        lambda^2
                end
                a = middle[e]
                gx[(j, a)] .+= gy * transpose(X[(a, k)])
                gx[(a, k)] .+= transpose(X[(j, a)]) * gy
            end
        end
        for (key, indices) in model.indices
            gradient[indices] .+= gx[key]
        end
    else
        @inbounds for (row, col, i, j, k, m, c) in model.quartic
            factor = c * l[row] * r[col] / lambda^2
            gradient[i] += factor * x[j] * x[k] * x[m]
            gradient[j] += factor * x[i] * x[k] * x[m]
            gradient[k] += factor * x[i] * x[j] * x[m]
            gradient[m] += factor * x[i] * x[j] * x[k]
        end
    end
    @inbounds for (row, col, i, j, c) in model.t.entries
        factor =
            c * (
                (l[row] * right_response[col] + left_response[row] * r[col]) / lambda^2 -
                2value / lambda * l[row] * r[col]
            )
        gradient[i] += factor * x[j]
        gradient[j] += factor * x[i]
    end
    return value, gradient
end

"Positive local tangent metric from the transfer fixed points, for search only."
function transfer_preconditioner(model, x; floor = 1e-4)
    W = ST.singlet_transfer(model.t, x)
    W = symmetric_restriction(W, model.symmetric_columns)
    lambda, r, l = fixed_points(W)
    r = symmetric_expansion(r, model.symmetric_columns, length(model.t.singlets))
    l = symmetric_expansion(l, model.symmetric_columns, length(model.t.singlets))
    sum(r) < 0 && (r = -r)
    l /= dot(l, r)
    function inverse_positive(A, reference)
        s = eigen(Symmetric((A+transpose(A))/2))
        # Newly added multiplets can have weight below eigensolver roundoff.
        # A common floor also prevents such blocks dominating a search step.
        return s.vectors*Diagonal(inv.(max.(s.values, floor*reference)))*transpose(s.vectors)
    end
    left, right = Dict{Int,Matrix{Float64}}(), Dict{Int,Matrix{Float64}}()
    for j in model.js
        d = length(model.groups[j])
        left[j] = transpose(reshape(l[model.blocks[j]], d, d))/(2j+1)
        right[j] = transpose(reshape(r[model.blocks[j]], d, d))
    end
    for blocks in (left, right)
        reference = maximum(opnorm(A) for A in values(blocks))
        reference > 0 || error("zero fixed-point metric")
        for j in model.js
            blocks[j] = inverse_positive(blocks[j], reference)
        end
    end
    function apply_metric(g)
        result = zero(g)
        for ((j, k), indices) in model.indices
            result[indices] = lambda*left[j]*g[indices]*right[k]
        end
        return result
    end
    return apply_metric
end

function optimize_lbfgs(
    model,
    initial;
    iterations = 2000,
    tolerance = 1e-8,
    memory = 20,
    precondition = false,
)
    x = initial / norm(initial)
    value, gradient = energy_gradient(model, x)
    inverse_metric = precondition ? transfer_preconditioner(model, x) : identity
    ss, ys, rhos = Vector{Float64}[], Vector{Float64}[], Float64[]
    reason = "iteration-limit"
    completed = 0
    for iteration in 1:iterations
        norm(gradient, Inf) < tolerance && (reason = "gradient-tolerance"; break)
        direction = copy(gradient)
        alpha = zeros(length(ss))
        for i in length(ss):-1:1
            alpha[i] = rhos[i] * dot(ss[i], direction)
            direction -= alpha[i] * ys[i]
        end
        direction = inverse_metric(direction)
        !isempty(ss) &&
            (direction *= dot(last(ss), last(ys)) / dot(last(ys), inverse_metric(last(ys))))
        for i in eachindex(ss)
            direction += ss[i] * (alpha[i] - rhos[i] * dot(ys[i], direction))
        end
        direction = -direction
        if dot(direction, gradient) >= 0
            empty!(ss)
            empty!(ys)
            empty!(rhos)
            direction = -inverse_metric(gradient)
        end
        norm(direction) > 1 && (direction /= norm(direction))
        step = 1.0
        candidate = x + step * direction
        fv, gn = energy_gradient(model, candidate)
        while (!isfinite(fv) || fv > value + 1e-4 * step * dot(gradient, direction)) && step > 1e-10
            step /= 2
            candidate = x + step * direction
            fv, gn = energy_gradient(model, candidate)
        end
        step <= 1e-10 && (reason = "line-search-floor"; break)
        s, y = candidate - x, gn - gradient
        sy = dot(s, y)
        if sy > 1e-16 * norm(s) * norm(y)
            push!(ss, s)
            push!(ys, y)
            push!(rhos, inv(sy))
            length(ss) > memory && (popfirst!(ss); popfirst!(ys); popfirst!(rhos))
        end
        x, value, gradient = candidate, fv, gn
        completed = iteration
    end
    x /= norm(x)
    value, gradient = energy_gradient(model, x)
    return (;
        x,
        energy = value,
        gradient_norm = norm(gradient, Inf),
        iterations = completed,
        reason,
    )
end

"BFGS with analytic gradients; canonicalization is applied between batches."
function optimize(model, initial; iterations = 1000, tolerance = 1e-9)
    x = initial / norm(initial)
    value, gradient = energy_gradient(model, x)
    inverse_hessian = Matrix{Float64}(I, length(x), length(x))
    completed = 0
    reason = "iteration-limit"
    for iteration in 1:iterations
        if norm(gradient, Inf) < tolerance
            reason = "gradient-tolerance"
            break
        end
        direction = -inverse_hessian * gradient
        if dot(direction, gradient) >= -1e-16 * norm(direction) * norm(gradient)
            direction = -gradient
            inverse_hessian = Matrix{Float64}(I, length(x), length(x))
        end
        norm(direction) > 1 && (direction /= norm(direction))
        step = 1.0
        candidate = x + step * direction
        fv, gn = energy_gradient(model, candidate)
        while (!isfinite(fv) || fv > value + 1e-4 * step * dot(gradient, direction)) && step > 1e-10
            step /= 2
            candidate = x + step * direction
            fv, gn = energy_gradient(model, candidate)
        end
        if step <= 1e-10
            reason = "line-search-floor"
            break
        end
        s, y = candidate - x, gn - gradient
        sy = dot(s, y)
        if sy > 1e-16 * norm(s) * norm(y)
            rho = inv(sy)
            Q = I - rho * s * transpose(y)
            inverse_hessian = Q * inverse_hessian * transpose(Q) + rho * s * transpose(s)
        end
        x, value, gradient = candidate, fv, gn
        completed = iteration
    end
    x /= norm(x)
    value, gradient = energy_gradient(model, x)
    return (;
        x,
        energy = value,
        gradient_norm = norm(gradient, Inf),
        iterations = completed,
        reason,
    )
end

function embed(old_t, old_x, new_t; noise = 0.005, seed = 1)
    rng = MersenneTwister(seed)
    x = noise * randn(rng, length(new_t.pairs))
    mapping = [
        findall(==(j), new_t.spins)[count(==(j), old_t.spins[1:a])] for
        (a, j) in enumerate(old_t.spins)
    ]
    lookup = Dict(pair => i for (i, pair) in enumerate(new_t.pairs))
    for (i, (a, b)) in enumerate(old_t.pairs)
        x[lookup[(mapping[a], mapping[b])]] = old_x[i]
    end
    return x / norm(x)
end

end
