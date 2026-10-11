"""
    ThermalMPOGap.CertifiedMPO.Norms

Rigorous operator and Hilbert--Schmidt norm witnesses for exact rational MPOs: local maps,
gauges, interval Cholesky acceptance, squared-tail bounds, and the sparse exact core encoding.
Moved verbatim from spin2 `Spin2MPONorms`; `encode_fraction`/`parse_fraction` are the
`ExactArithmetic` ones (identical code), re-bound here so `Norms.encode_fraction` still works.
"""
module Norms

using LinearAlgebra
import Nemo
using ..CertifiedMPO: F, nemo_array, encode_fraction, parse_fraction

"Numerical candidate generation with a QR fallback for a failed divide-and-conquer SVD."
function numerical_svd(A; full = false)
    @assert all(isfinite, A)
    try
        return svd(A; full)
    catch err
        err isa LinearAlgebra.LAPACKException || rethrow()
        return svd(A; alg = LinearAlgebra.QRIteration(), full)
    end
end

"Dyadic upper square root, accepted by an integer square comparison."
function upper_sqrt(x::F; bits = 80)
    @assert x >= 0 && bits >= 0
    denominator = big(2)^bits
    square = cld(numerator(x)*denominator^2, Base.denominator(x))
    root = isqrt(square)
    root^2 < square && (root += 1)
    result = root//denominator
    @assert result^2 >= x
    return result
end

"Scale before taking a dyadic upper square root, preserving relative accuracy for tiny residuals."
function relative_upper_sqrt(x::F; bits = 80)
    @assert x >= 0
    iszero(x) && return F(0)
    exponent = fld(ndigits(numerator(x); base = 2)-ndigits(denominator(x); base = 2), 2)
    scale = exponent >= 0 ? F(big(2)^exponent) : F(1, big(2)^(-exponent))
    return scale*upper_sqrt(x/scale^2; bits)
end

"Round a nonnegative rational upward to a dyadic with bounded relative precision."
function relative_rational_upper(x::F; bits = 128)
    @assert x >= 0 && bits >= 1
    iszero(x) && return F(0)
    exponent = ndigits(numerator(x); base = 2)-ndigits(denominator(x); base = 2)-bits
    scale = exponent >= 0 ? F(big(2)^exponent) : F(1, big(2)^(-exponent))
    scaled = x/scale
    result = cld(numerator(scaled), denominator(scaled))*scale
    @assert x <= result < x+scale
    return result
end

"Quantize a real MPO candidate; this defines its exact rational coefficients."
function quantize(cores; bits = 48)
    denominator = big(2)^bits
    @assert all(all(isreal, core) for core in cores)
    return [[round(BigInt, real(x)*denominator)//denominator for x in core] for core in cores]
end

"Bound one partial-operator family through a core (left, output, input, right)."
function propagate(core::Array{F,4}, previous; bits = 80)
    l, dout, din, r = size(core)
    @assert l == length(previous) && all(>=(0), previous)
    result = zeros(F, r)
    for a in 1:r
        envelope = zeros(F, dout, din)
        for s in 1:dout, t in 1:din, b in 1:l
            envelope[s, t] += abs(core[b, s, t, a])*previous[b]
        end
        # If ||X_st|| <= C_st, then ||[X_st]|| <= ||C||_2.
        one_norm = maximum(sum(envelope; dims = 1))
        infinity_norm = maximum(sum(envelope; dims = 2))
        result[a] = upper_sqrt(one_norm*infinity_norm; bits)
    end
    return result
end

"Rigorous operator-norm witnesses from both directions and every MPO cut."
function witness(cores; bits = 80)
    @assert !isempty(cores) && size(first(cores), 1) == size(last(cores), 4) == 1
    @assert all(size(cores[j], 4) == size(cores[j+1], 1) for j in 1:(length(cores)-1))
    left = Vector{F}[F[1]]
    for core in cores
        push!(left, propagate(core, last(left); bits))
    end
    right = [F[] for _ in 1:(length(cores)+1)]
    right[end] = F[1]
    for j in length(cores):-1:1
        right[j] = propagate(permutedims(cores[j], (4, 2, 3, 1)), right[j+1]; bits)
    end
    cuts = [sum(left[j] .* right[j]) for j in eachindex(left)]
    return (; upper = minimum(cuts), cuts, left, right)
end

"View a core as an operator from (left,input) to (right,output)."
function local_operator(core)
    l, dout, din, r = size(core)
    return reshape(permutedims(core, (2, 4, 3, 1)), dout*r, din*l)
end

"Contract two adjacent cores exactly, preserving their physical operator."
function merge_cores(first, second; field = Nemo.QQ)
    l, do1, di1, mid = size(first)
    @assert size(second, 1) == mid
    _, do2, di2, r = size(second)
    A = Nemo.matrix(field, reshape(first, l*do1*di1, mid))
    B = Nemo.matrix(field, reshape(second, mid, do2*di2*r))
    product = A*B
    raw = reshape(nemo_array(product, field), l, do1, di1, do2, di2, r)
    return reshape(permutedims(raw, (1, 4, 2, 5, 3, 6)), l, do1*do2, di1*di2, r)
end

"Exact positive rational gauges, proposed by a numerical logarithmic search."
function rational_gauges(cores, x; bits = 40)
    @assert all(isfinite, x) && bits >= 0
    dimensions = [size(c, 4) for c in cores[1:(end-1)]]
    @assert length(x) == sum(dimensions)
    denominator = big(2)^bits
    result = Vector{F}[F[1]]
    offset = 0
    for d in dimensions
        z = x[(offset+1):(offset+d)]
        # A common factor on one bond cancels in the product of local norms.
        weights = exp.(z .- maximum(z))
        push!(result, [max(big(1), round(BigInt, w*denominator))//denominator for w in weights])
        offset += d
    end
    push!(result, F[1])
    return result
end

"Apply rational vector/diagonal or full-matrix gauges with exact cancellation."
function gauge_transform(cores, gauges; interval_precision = nothing)
    @assert length(gauges) == length(cores)+1
    if first(gauges) isa AbstractMatrix
        @assert all(g -> g isa AbstractMatrix, gauges)
        @assert first(gauges) == last(gauges) == ones(F, 1, 1)
        matrices = [Nemo.matrix(Nemo.QQ, g) for g in gauges]
        inverses = inv.(matrices)
        for (g, inverse) in zip(matrices, inverses)
            @assert g*inverse == Nemo.identity_matrix(Nemo.QQ, Nemo.nrows(g))
        end
        RR = isnothing(interval_precision) ? nothing : Nemo.ArbField(interval_precision)
        interval_matrices = isnothing(RR) ? nothing : [Nemo.matrix(RR, g) for g in gauges]
        interval_inverses =
            isnothing(RR) ? nothing : [Nemo.matrix(RR, nemo_array(g)) for g in inverses]
        result = isnothing(RR) ? Array{F,4}[] : Any[]
        for (j, core) in enumerate(cores)
            l, d, di, r = size(core)
            @assert size(gauges[j]) == (l, l) && size(gauges[j+1]) == (r, r)
            # The terminal core remains exact for the final bridge equality.
            # Elsewhere only an enclosure is needed for residual/norm checks.
            field = isnothing(RR) || j == length(cores) ? Nemo.QQ : RR
            g = field === Nemo.QQ ? matrices[j] : interval_matrices[j]
            inverse = field === Nemo.QQ ? inverses[j+1] : interval_inverses[j+1]
            left = g*Nemo.matrix(field, reshape(core, l, d*di*r))
            right = Nemo.matrix(field, reshape(nemo_array(left, field), l*d*di, r))*inverse
            push!(result, reshape(nemo_array(right, field), l, d, di, r))
        end
        return result
    end
    @assert all(g -> g isa AbstractVector, gauges)
    @assert first(gauges) == last(gauges) == F[1]
    @assert all(all(>(0), g) for g in gauges)
    result = Array{F,4}[]
    for (j, core) in enumerate(cores)
        l, dout, din, r = size(core)
        @assert length(gauges[j]) == l && length(gauges[j+1]) == r
        push!(
            result,
            [
                core[b, s, t, a]*gauges[j][b]/gauges[j+1][a] for
                b in 1:l, s in 1:dout, t in 1:din, a in 1:r
            ],
        )
    end
    return result
end

"One exact transformed core for an inconclusive interval-residual fallback."
function exact_gauge_core(core, left_gauge, right_gauge)
    l, d, di, r = size(core)
    if left_gauge isa AbstractVector
        @assert length(left_gauge) == l && length(right_gauge) == r
        return core .* reshape(left_gauge, l, 1, 1, 1) ./ reshape(right_gauge, 1, 1, 1, r)
    end
    @assert size(left_gauge) == (l, l) && size(right_gauge) == (r, r)
    left_matrix, right_matrix = Nemo.matrix(Nemo.QQ, left_gauge), Nemo.matrix(Nemo.QQ, right_gauge)
    inverse = inv(right_matrix)
    @assert right_matrix*inverse == Nemo.identity_matrix(Nemo.QQ, r)
    left = left_matrix*Nemo.matrix(Nemo.QQ, reshape(core, l, d*di*r))
    right = Nemo.matrix(Nemo.QQ, reshape(nemo_array(left), l*d*di, r))*inverse
    return reshape(nemo_array(right), l, d, di, r)
end

"Portable gauges, retaining the original vector encoding for diagonal witnesses."
encode_gauges(gauges) = [
    g isa AbstractMatrix ?
    Dict("dimensions"=>collect(size(g)), "entries"=>encode_fraction.(vec(g))) :
    encode_fraction.(g) for g in gauges
]

decode_gauges(encoded) = [
    g isa AbstractVector ? parse_fraction.(g) :
    reshape(parse_fraction.(g["entries"]), Int.(g["dimensions"])...) for g in encoded
]

"Build an outward-rounded Gram matrix in the smaller local space."
function local_gram(W::AbstractMatrix; precision = 192)
    RR = Nemo.ArbField(precision)
    M = Nemo.matrix(RR, W)
    return size(W, 1) <= size(W, 2) ? M*transpose(M) : transpose(M)*M
end

"Accept a local spectral norm upper bound only after rigorous Arb Cholesky."
function verify_gram_bound(gram, alpha::F; metric = nothing)
    alpha > 0 || return false
    RR = Nemo.base_ring(gram)
    C = -gram
    for i in 1:Nemo.nrows(C)
        C[i, i] += RR(alpha)^2*(isnothing(metric) ? RR(1) : RR(metric[i]))
    end
    try
        Nemo.cholesky(C)
        return true
    catch err
        # arb_mat_cho returns failure if positive definiteness is not certain.
        if err isa ErrorException && err.msg == "Matrix is not positive definite"
            return false
        end
        rethrow()
    end
end

"Rational weighted Gram pair for the norm from column metric to row metric."
function weighted_gram(W; row_metric = nothing, col_metric = nothing, precision = 192)
    if isnothing(row_metric) && isnothing(col_metric)
        return local_gram(W; precision), nothing
    end
    rows = isnothing(row_metric) ? ones(F, size(W, 1)) : F.(row_metric)
    cols = isnothing(col_metric) ? ones(F, size(W, 2)) : F.(col_metric)
    @assert length(rows) == size(W, 1) && length(cols) == size(W, 2)
    @assert all(>(0), rows) && all(>(0), cols)
    RR = Nemo.ArbField(precision)
    M = Nemo.matrix(RR, W)
    scaled = deepcopy(M)
    if size(W, 1) <= size(W, 2)
        for j in eachindex(cols), i in eachindex(rows)
            scaled[i, j] /= RR(cols[j])
        end
        return scaled*transpose(M), inv.(rows)
    end
    for j in eachindex(cols), i in eachindex(rows)
        scaled[i, j] *= RR(rows[i])
    end
    return transpose(M)*scaled, cols
end

function verify_local_bound(
    W::AbstractMatrix,
    alpha::F;
    precision = 192,
    row_metric = nothing,
    col_metric = nothing,
)
    alpha == 0 && return all(iszero, W)
    alpha < 0 && return false
    gram, metric = weighted_gram(W; precision, row_metric, col_metric)
    return verify_gram_bound(gram, alpha; metric)
end

"Numerical proposals are independently accepted or rejected by ball arithmetic."
function local_norm_upper(
    W::AbstractMatrix;
    precision = 192,
    bits = 48,
    row_metric = nothing,
    col_metric = nothing,
)
    all(iszero, W) && return F(0)
    candidate = Float64.(W)
    if !isnothing(row_metric)
        candidate .*= sqrt.(Float64.(row_metric))
    end
    if !isnothing(col_metric)
        candidate ./= transpose(sqrt.(Float64.(col_metric)))
    end
    estimate = try
        opnorm(candidate)
    catch err
        err isa LinearAlgebra.LAPACKException || rethrow()
        first(svd(candidate; alg = LinearAlgebra.QRIteration()).S)
    end
    @assert isfinite(estimate) && estimate > 0 && bits >= 0
    gram, metric = weighted_gram(W; precision, row_metric, col_metric)
    denominator = big(2)^bits
    exponent = last(frexp(estimate))-1
    scale = exponent >= 0 ? F(big(2)^exponent) : F(1, big(2)^(-exponent))
    for attempt in 0:12
        proposal = estimate*(1+10.0^(attempt-10))
        alpha = scale*((ceil(BigInt, (proposal/Float64(scale))*denominator)+1)//denominator)
        verify_gram_bound(gram, alpha; metric) && return alpha
    end
    error("No local norm bound passed the interval Cholesky check")
end

"Physical Hilbert--Schmidt row metric, in (left,output,input) order."
function hs_row_metric(left, d, weights)
    isnothing(weights) && return ones(F, left*d*d)
    @assert length(weights) == d && all(>(0), weights)
    return vec([F(weights[t])/F(weights[s]) for b in 1:left, s in 1:d, t in 1:d])
end

"A local rational map and its physical input/output metrics."
function core_map(core, kind, weights = nothing)
    l, d, di, r = size(core)
    @assert d == di
    @assert isnothing(weights) || (length(weights)==d && all(>(0), weights))
    if kind == :operator
        W = local_operator(core)
        rows = isnothing(weights) ? nothing : repeat(inv.(F.(weights)), r)
        cols = isnothing(weights) ? nothing : repeat(inv.(F.(weights)), l)
    elseif kind == :left_hs
        W = left_hs_map(core)
        rows = isnothing(weights) ? nothing : hs_row_metric(l, d, weights)
        cols = nothing
    elseif kind == :right_hs
        W = right_hs_map(core)
        rows = isnothing(weights) ? nothing : repeat(hs_row_metric(1, d, weights), r)
        cols = nothing
    else
        error("Unknown local map kind")
    end
    return W, rows, cols
end

"Physical adjoint in the rational polynomial basis; virtual order is unchanged."
function adjoint_cores(cores; weights = nothing)
    return [
        [
            (isnothing(weights) ? F(1) : F(weights[s])/F(weights[t]))*c[b, t, s, a] for
            b in axes(c, 1), s in axes(c, 3), t in axes(c, 2), a in axes(c, 4)
        ] for c in cores
    ]
end

"Exact spatial reflection. This is unitary on the full physical Hilbert space."
reflect_cores(cores) = [permutedims(c, (4, 2, 3, 1)) for c in reverse(cores)]

function core_norm_upper(core, kind; weights = nothing, kwargs...)
    W, rows, cols = core_map(core, kind, weights)
    return local_norm_upper(W; row_metric = rows, col_metric = cols, kwargs...)
end

function verify_core_bound(core, kind, alpha; weights = nothing, kwargs...)
    W, rows, cols = core_map(core, kind, weights)
    return verify_local_bound(W, alpha; row_metric = rows, col_metric = cols, kwargs...)
end

"Global MPO norm bound from a sequence of rigorously bounded local maps."
function local_product_witness(cores, gauges; precision = 192, bits = 48, weights = nothing)
    transformed = gauge_transform(cores, gauges)
    bounds = [core_norm_upper(c, :operator; precision, bits, weights) for c in transformed]
    return (; upper = prod(bounds), bounds)
end

"Independent replay of proposed local bounds, without an SVD or optimizer."
function verify_local_product(
    cores,
    gauges,
    bounds;
    precision = 192,
    weights = nothing,
    interval_precision = nothing,
)
    length(bounds) == length(cores) || return false
    transformed = gauge_transform(cores, gauges; interval_precision)
    return all(
        verify_core_bound(c, :operator, alpha; precision, weights) for
        (c, alpha) in zip(transformed, bounds)
    )
end

"Sparse exact core encoding, independent of floating-point JSON semantics."
function encode_cores(cores)
    return [
        Dict(
            "dimensions"=>collect(size(c)),
            "entries"=>[[i, encode_fraction(c[i])] for i in eachindex(c) if !iszero(c[i])],
        ) for c in cores
    ]
end

function decode_cores(encoded)
    result = Array{F,4}[]
    for row in encoded
        dimensions = Int.(row["dimensions"])
        @assert length(dimensions) == 4 && all(>(0), dimensions)
        c = zeros(F, dimensions...)
        seen = Set{Int}()
        for pair in row["entries"]
            i = Int(pair[1])
            @assert 1 <= i <= length(c) && !(i in seen)
            push!(seen, i)
            c[i] = parse_fraction(pair[2])
        end
        push!(result, c)
    end
    return result
end

"A right-tail coefficient map to Hilbert--Schmidt operators."
function right_hs_map(core)
    l, dout, din, r = size(core)
    return reshape(permutedims(core, (2, 3, 4, 1)), dout*din*r, l)
end

"A left-prefix coefficient map to Hilbert--Schmidt operators."
left_hs_map(core) = reshape(core, prod(size(core)[1:3]), size(core, 4))

"Bounds on the coefficient maps of all squared MPO tails, without a four-layer contraction."
function squared_tail_witness(
    cores;
    gauges = nothing,
    precision = 192,
    cap_frobenius = true,
    blocksize = 1,
    weights = nothing,
    adjoint_first = false,
    return_operator_cores = false,
    arithmetic = :rational,
)
    @assert blocksize in (1, 2)
    @assert arithmetic in (:rational, :interval)
    first_copy = adjoint_first ? adjoint_cores(cores; weights) : cores
    operator_cores =
        isnothing(gauges) ? first_copy :
        gauge_transform(
            first_copy,
            gauges;
            interval_precision = arithmetic == :interval ? precision : nothing,
        )
    local_operator_bounds =
        [core_norm_upper(c, :operator; weights, precision) for c in operator_cores]
    local_hs_bounds = [core_norm_upper(c, :right_hs; weights, precision) for c in cores]
    block_bounds =
        blocksize == 2 ?
        [
            core_norm_upper(
                merge_cores(
                    operator_cores[j],
                    operator_cores[j+1];
                    field = arithmetic == :interval ? Nemo.ArbField(precision) : Nemo.QQ,
                ),
                :operator;
                weights = isnothing(weights) ? nothing : kron(weights, weights),
                precision,
            ) for j in 1:(length(cores)-1)
        ] : F[]
    operator_tail = ones(F, length(cores)+1)
    hs_tail = ones(F, length(cores)+1)
    for j in length(cores):-1:1
        operator_tail[j] = operator_tail[j+1]*local_operator_bounds[j]
        if blocksize == 2 && j < length(cores)
            operator_tail[j] = min(operator_tail[j], block_bounds[j]*operator_tail[j+2])
        end
        hs_tail[j] = hs_tail[j+1]*local_hs_bounds[j]
        if cap_frobenius
            # The first-copy tails are g_a R_a. The HS coefficient-map
            # bound gives sum_a ||g_a R_a||_F^2 <= gamma^2 sum_a g_a^2.
            weight2 = isnothing(gauges) ? F(size(cores[j], 1)) : sum(abs2, gauges[j])
            operator_tail[j] = min(operator_tail[j], hs_tail[j]*relative_upper_sqrt(weight2))
        end
    end
    # For R(v)=sum_a v_a R_a and K(V)=sum_ab V_ab R_a R_b,
    # K(V)=[R_1 ... R_D] col_a(sum_b V_ab R_b).
    # Therefore ||K||_(F <- F) <= ||[R_1 ... R_D]||_op ||R||_(F <- 2).
    certificate = (;
        upper = operator_tail .* hs_tail,
        operator_tail,
        hs_tail,
        local_operator_bounds,
        local_hs_bounds,
        frobenius_capped = cap_frobenius,
        block_bounds,
    )
    return return_operator_cores ? (; certificate, operator_cores) : certificate
end

end
