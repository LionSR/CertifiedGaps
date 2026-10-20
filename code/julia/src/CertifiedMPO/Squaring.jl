"""
    CertifiedGaps.CertifiedMPO.Squaring

Certified left-to-right compression of an MPO square B² (or reflected Gram B'B), its exact
replay, portable witnesses, and Gram-flow error propagation. Moved verbatim from spin2
`Spin2MPOSquaring`; `MN` is `CertifiedMPO.Norms` instead of a private include copy.
"""
module Squaring

using LinearAlgebra
import Nemo
using ..CertifiedMPO: nemo_array
import ..Norms
const MN = Norms
const F = MN.F
import ..ChargeBlocks
const CB = ChargeBlocks

function matrix_product(A, B, field)
    C = Nemo.matrix(field, A)*Nemo.matrix(field, B)
    return nemo_array(C, field)
end

rational_product(A, B) = matrix_product(A, B, Nemo.QQ)

"Contract a bridge with a product core, allowing unequal auxiliary dimensions."
function product_bridge(bridge, core, second; field = Nemo.QQ)
    l, d, mid, r = size(core)
    lb, mid2, di, rb = size(second)
    @assert mid == mid2 && size(bridge, 2) == l && size(bridge, 3) == lb
    chi = size(bridge, 1)
    X = reshape(permutedims(bridge, (1, 3, 2)), chi*lb, l)
    Y = reshape(core, l, d*mid*r)
    A = reshape(matrix_product(X, Y, field), chi, lb, d, mid, r)
    X2 = reshape(permutedims(A, (1, 3, 5, 2, 4)), chi*d*r, lb*mid)
    Y2 = reshape(second, lb*mid, di*rb)
    B = reshape(matrix_product(X2, Y2, field), chi, d, r, di, rb)
    return permutedims(B, (1, 2, 4, 3, 5))
end

"Contract an incoming bridge with one squared MPO core, without constructing that core."
square_bridge(bridge, core, second = core; field = Nemo.QQ) =
    product_bridge(bridge, core, second; field)

"Dyadic rounding with a scale chosen from the candidate; no accuracy is assumed."
function relative_quantize(A; bits = 40)
    largest = maximum(abs, A)
    iszero(largest) && return zeros(F, size(A))
    @assert isfinite(largest)
    exponent = exponent2(largest)
    scale = exponent >= 0 ? F(big(2)^exponent) : F(1, big(2)^(-exponent))
    denominator = big(2)^bits
    return [round(BigInt, (x/Float64(scale))*denominator)//denominator*scale for x in A]
end

exponent2(x) = last(frexp(x))-1

"Upper Frobenius norm of an exact matrix or an outward-rounded enclosure."
function frobenius_upper(A, field; metric = nothing)
    weights = isnothing(metric) ? nothing : field === Nemo.QQ ? metric : field.(metric)
    squared = isnothing(weights) ? sum(abs2, A) : sum(abs2.(A) .* weights)
    if field !== Nemo.QQ
        endpoint = BigFloat(squared, RoundUp)
        @assert isfinite(endpoint) && endpoint >= 0
        squared = F(endpoint)
    end
    return MN.relative_upper_sqrt(F(squared))
end

"Backward error recurrence retaining certified projection cross terms."
function projection_error_upper(rows, tail_upper; bits = 128)
    @assert length(tail_upper) == length(rows)+2
    error = F(0)
    for j in reverse(eachindex(rows))
        row = rows[j]
        delta, alpha, theta = row.residual_bound, row.local_prefix_bound, row.projection_bound
        @assert min(delta, alpha, theta, tail_upper[j+1]) >= 0
        # ||U'*E||_F <= ||U|| * ||E||_F is a second, independent bound.
        theta = min(theta, alpha*delta)
        local_error = delta*tail_upper[j+1]
        squared = local_error^2+(alpha*error)^2+2theta*tail_upper[j+1]*error
        rounded = MN.relative_rational_upper(MN.relative_upper_sqrt(squared); bits)
        error = min(rounded, local_error+alpha*error)
    end
    return error
end

"Change only auxiliary scales, keeping every compressed output core exact."
function rescale_square_gauges(certificate, scales)
    @assert !isnothing(certificate.gauges)
    @assert length(scales) == length(certificate.gauges) == length(certificate.cores)+1
    scales = [F.(s) for s in scales]
    @assert first(scales) == last(scales) == F[1]
    @assert all(s -> all(>(0), s), scales)
    gauges = map(certificate.gauges, scales) do g, s
        @assert length(s) == size(g, 1)
        s .* g
    end
    bridges = [B ./ reshape(scales[j+1], 1, :, 1) for (j, B) in enumerate(certificate.bridges)]
    # G_new = D*G_old, hence R_new = R_old*D^-1 on the
    # first input's auxiliary leg. Local products cancel D exactly.
    # New residual and tail bounds must be generated and verified separately.
    return (; cores = certificate.cores, bridges, gauges)
end

"One exact residual certificate for a left-to-right SVD compression of B squared."
function square_candidate(
    cores;
    maxdim = 8,
    bits = 40,
    precision = 192,
    gauges = nothing,
    blocksize = 1,
    weights = nothing,
    adjoint_first = false,
    candidate_arithmetic = :rational,
    error_method = :triangle,
    project_bridge = false,
    fixed_compression = nothing,
    physical_charges = nothing,
    progress = nothing,
)
    @assert !isempty(cores) && size(first(cores), 1) == size(last(cores), 4) == 1
    @assert candidate_arithmetic in (:rational, :interval)
    @assert error_method in (:triangle, :projection)
    if !isnothing(fixed_compression)
        @assert length(fixed_compression.cores) == length(cores)
        @assert length(fixed_compression.bridges) == length(cores)-1
    end
    RR = Nemo.ArbField(precision)
    isnothing(progress) || progress((; stage = :tail_bounds, site = 0, total = length(cores)))
    prepared = MN.squared_tail_witness(
        cores;
        gauges,
        precision,
        blocksize,
        weights,
        adjoint_first,
        return_operator_cores = true,
        arithmetic = candidate_arithmetic,
    )
    tail, operator_cores = prepared.certificate, prepared.operator_cores
    charge_a =
        isnothing(physical_charges) ? nothing : CB.infer_charges(operator_cores, physical_charges)
    charge_b = isnothing(physical_charges) ? nothing : CB.infer_charges(cores, physical_charges)
    prefix_charges = [0]
    @assert isnothing(physical_charges) || isnothing(fixed_compression)
    isnothing(progress) || progress((; stage = :compression, site = 0, total = length(cores)))
    bridge = ones(F, 1, 1, 1)
    result = Array{F,4}[]
    bridges = Array{F,3}[]
    prefix = F(1)
    error = F(0)
    rows = Any[]
    for (j, core) in enumerate(cores)
        # Keep the terminal core exact. Interior bridge products may be
        # enclosed instead: the SVD proposes exact dyadics, and the residual
        # bound includes the complete interval enclosure.
        field = candidate_arithmetic == :interval && j < length(cores) ? RR : Nemo.QQ
        Z =
            isnothing(physical_charges) ? square_bridge(bridge, operator_cores[j], core; field) :
            CB.product_bridge(
                bridge,
                operator_cores[j],
                core,
                prefix_charges,
                charge_a[j],
                charge_a[j+1],
                charge_b[j],
                charge_b[j+1],
                physical_charges;
                field,
            )
        chi, d, _, r, _ = size(Z)
        M = reshape(Z, chi*d*d, r*r)
        if j == length(cores)
            isnothing(fixed_compression) ||
                @assert reshape(fixed_compression.cores[j], size(M)) == M
            push!(result, reshape(M, chi, d, d, 1))
            isnothing(progress) ||
                progress((; stage = :compression, site = j, total = length(cores)))
            break
        end
        row_metric = MN.hs_row_metric(chi, d, weights)
        root_metric = sqrt.(Float64.(row_metric))
        U, R = if isnothing(fixed_compression)
            labels =
                isnothing(physical_charges) ? nothing :
                CB.product_labels(prefix_charges, charge_a[j+1], charge_b[j+1], physical_charges)
            decomposition =
                isnothing(labels) ? MN.numerical_svd(Float64.(M) .* root_metric) :
                CB.truncated_svd(
                    Float64.(M) .* root_metric,
                    labels...;
                    maxdim,
                    svd_function = MN.numerical_svd,
                )
            k = min(maxdim, length(decomposition.S))
            proposed_left = relative_quantize(decomposition.U[:, 1:k] ./ root_metric; bits)
            # Fit against the actual rounded columns, not their unrounded SVD
            # predecessors. The residual and projection defect still certify
            # every arithmetic error in this purely numerical proposal.
            proposed_right =
                project_bridge ?
                (
                    isnothing(labels) ?
                    (Float64.(proposed_left) .* root_metric) \ (Float64.(M) .* root_metric) :
                    CB.fit_bridge(
                        Float64.(proposed_left) .* root_metric,
                        Float64.(M) .* root_metric,
                        labels[1],
                        decomposition.charges,
                        labels[2],
                    )
                ) : Diagonal(decomposition.S[1:k])*decomposition.Vt[1:k, :]
            isnothing(labels) || (prefix_charges = decomposition.charges)
            proposed_left, relative_quantize(proposed_right; bits)
        else
            k = size(fixed_compression.cores[j], 4)
            reshape(fixed_compression.cores[j], chi*d*d, k),
            reshape(fixed_compression.bridges[j], k, r*r)
        end
        k = size(U, 2)
        residual = M-matrix_product(U, R, field)
        residual_bound = frobenius_upper(residual, field; metric = row_metric)
        projection_bound =
            error_method == :projection ?
            frobenius_upper(matrix_product(transpose(U .* row_metric), residual, field), field) :
            F(0)
        term = prefix*residual_bound*tail.upper[j+1]
        error += term
        next = reshape(U, chi, d, d, k)
        push!(result, next)
        local_prefix_bound = MN.core_norm_upper(next, :left_hs; weights, precision)
        push!(
            rows,
            (;
                site = j,
                residual_bound,
                prefix_bound = prefix,
                squared_tail_bound = tail.upper[j+1],
                error_term = term,
                local_prefix_bound,
                projection_bound,
            ),
        )
        prefix *= local_prefix_bound
        bridge = reshape(R, k, r, r)
        push!(bridges, bridge)
        isnothing(progress) || progress((; stage = :compression, site = j, total = length(cores)))
    end
    error = error_method == :projection ? projection_error_upper(rows, tail.upper) : error
    return (;
        cores = result,
        bridges,
        error,
        rows,
        tail,
        gauges,
        weights,
        adjoint_first,
        error_method,
    )
end

"HS norm bound from already verified left-prefix maps and the small terminal core."
function output_hs_upper(certificate)
    C = last(certificate.cores)
    chi, d, di, r = size(C)
    @assert d == di && r == 1
    metric = MN.hs_row_metric(chi, d, certificate.weights)
    terminal2 = sum(abs2.(reshape(C, chi*d*d, 1)) .* metric)
    terminal = MN.relative_upper_sqrt(terminal2)
    return terminal*prod(row.local_prefix_bound for row in certificate.rows; init = F(1))
end

"Replay every local residual and interval bound without trusting the SVD search."
function verify_square_certificate(
    input,
    certificate;
    precision = 192,
    details = false,
    residual_arithmetic = :rational,
    progress = nothing,
)
    @assert residual_arithmetic in (:rational, :interval)
    RR = Nemo.ArbField(precision)
    n = length(input)
    @assert length(certificate.cores) == n
    @assert length(certificate.bridges) == length(certificate.rows) == n-1
    error_method = hasproperty(certificate, :error_method) ? certificate.error_method : :triangle
    @assert error_method in (:triangle, :projection)
    tail = certificate.tail
    gauges = hasproperty(certificate, :gauges) ? certificate.gauges : nothing
    weights = hasproperty(certificate, :weights) ? certificate.weights : nothing
    adjoint_first = hasproperty(certificate, :adjoint_first) && certificate.adjoint_first
    first_copy = adjoint_first ? MN.adjoint_cores(input; weights) : input
    operator_cores =
        isnothing(gauges) ? first_copy :
        MN.gauge_transform(
            first_copy,
            gauges;
            interval_precision = residual_arithmetic == :interval ? precision : nothing,
        )
    op_tail, hs_tail = ones(F, n+1), ones(F, n+1)
    for j in n:-1:1
        @assert MN.verify_core_bound(
            operator_cores[j],
            :operator,
            tail.local_operator_bounds[j];
            weights,
            precision,
        )
        @assert MN.verify_core_bound(
            input[j],
            :right_hs,
            tail.local_hs_bounds[j];
            weights,
            precision,
        )
        op_tail[j] = op_tail[j+1]*tail.local_operator_bounds[j]
        if hasproperty(tail, :block_bounds) && !isempty(tail.block_bounds) && j < n
            W = MN.merge_cores(
                operator_cores[j],
                operator_cores[j+1];
                field = residual_arithmetic == :interval ? RR : Nemo.QQ,
            )
            @assert MN.verify_core_bound(
                W,
                :operator,
                tail.block_bounds[j];
                precision,
                weights = isnothing(weights) ? nothing : kron(weights, weights),
            )
            op_tail[j] = min(op_tail[j], tail.block_bounds[j]*op_tail[j+2])
        end
        hs_tail[j] = hs_tail[j+1]*tail.local_hs_bounds[j]
        if hasproperty(tail, :frobenius_capped) && tail.frobenius_capped
            weight2 = isnothing(gauges) ? F(size(input[j], 1)) : sum(abs2, gauges[j])
            op_tail[j] = min(op_tail[j], hs_tail[j]*MN.relative_upper_sqrt(weight2))
        end
        isnothing(progress) || progress((; stage = :tail_bounds, site = n+1-j, total = n))
    end
    bridge = ones(F, 1, 1, 1)
    prefix, error = F(1), F(0)
    for j in 1:n
        field = residual_arithmetic == :interval && j < n ? RR : Nemo.QQ
        Z = square_bridge(bridge, operator_cores[j], input[j]; field)
        chi, d, _, r, _ = size(Z)
        M = reshape(Z, chi*d*d, r*r)
        C = certificate.cores[j]
        @assert size(C)[1:3] == (chi, d, d)
        if j == n
            @assert size(C, 4) == 1
            @assert reshape(C, size(M)) == M
            break
        end
        k = size(C, 4)
        next = certificate.bridges[j]
        @assert size(next) == (k, r, r)
        U, R = reshape(C, chi*d*d, k), reshape(next, k, r*r)
        row = certificate.rows[j]
        metric = MN.hs_row_metric(chi, d, weights)
        @assert row.residual_bound >= 0
        if residual_arithmetic == :interval
            residual = M-matrix_product(U, R, RR)
            squared = sum(abs2.(residual) .* RR.(metric))
            if !(RR(row.residual_bound)^2 > squared)
                # Zero residuals and bounds with exceptionally small slack may
                # need exact arithmetic; an inconclusive interval is not proof.
                exact_core =
                    isnothing(gauges) ? first_copy[j] :
                    MN.exact_gauge_core(first_copy[j], gauges[j], gauges[j+1])
                exact = reshape(square_bridge(bridge, exact_core, input[j]), size(M))
                rational_residual = exact-rational_product(U, R)
                @assert row.residual_bound^2 >= sum(abs2.(rational_residual) .* metric)
            end
        else
            residual = M-rational_product(U, R)
            @assert row.residual_bound^2 >= sum(abs2.(residual) .* metric)
        end
        if error_method == :projection
            @assert row.projection_bound >= 0
            projected = matrix_product(transpose(U .* metric), residual, field)
            if residual_arithmetic == :interval
                if !(RR(row.projection_bound)^2 > sum(abs2, projected))
                    exact_core =
                        isnothing(gauges) ? first_copy[j] :
                        MN.exact_gauge_core(first_copy[j], gauges[j], gauges[j+1])
                    exact = reshape(square_bridge(bridge, exact_core, input[j]), size(M))
                    rational_residual = exact-rational_product(U, R)
                    projected_exact = rational_product(transpose(U .* metric), rational_residual)
                    @assert row.projection_bound^2 >= sum(abs2, projected_exact)
                end
            else
                @assert row.projection_bound^2 >= sum(abs2, projected)
            end
        end
        @assert MN.verify_core_bound(C, :left_hs, row.local_prefix_bound; weights, precision)
        error += prefix*row.residual_bound*op_tail[j+1]*hs_tail[j+1]
        prefix *= row.local_prefix_bound
        bridge = next
        isnothing(progress) || progress((; stage = :residuals, site = j, total = n))
    end
    triangle_error = error
    error =
        error_method == :projection ? projection_error_upper(certificate.rows, op_tail .* hs_tail) :
        error
    @assert error == certificate.error
    # The gauges that bound local compression residuals need not give the
    # tightest global operator norm. Replay a separate product bound when
    # supplied, without changing any residual tail or its certified error.
    operator_norm_upper = op_tail[1]
    if hasproperty(certificate, :operator_certificate) &&
       !isnothing(certificate.operator_certificate)
        independent = certificate.operator_certificate
        @assert MN.verify_local_product(
            first_copy,
            independent.gauges,
            independent.bounds;
            precision,
            weights,
            interval_precision = residual_arithmetic == :interval ? precision : nothing,
        )
        operator_norm_upper = min(operator_norm_upper, prod(independent.bounds))
    end
    return details ?
           (;
        valid = true,
        operator_norm_upper,
        hs_norm_upper = hs_tail[1],
        output_hs_norm_upper = output_hs_upper(certificate),
        triangle_error_upper = triangle_error,
    ) : true
end

"Outward-rounded physical Hilbert--Schmidt inner product of two exact rational MPOs."
function hs_inner_ball(first, second; weights = nothing, precision = 192)
    @assert length(first) == length(second)
    RR = Nemo.ArbField(precision)
    G = Nemo.matrix(RR, ones(F, 1, 1))
    cache = Dict{Tuple{UInt,Int,Int},Any}()
    block(c, s, t) = get!(cache, (objectid(c), s, t)) do
        any(!iszero, c[:, s, t, :]) ? Nemo.matrix(RR, c[:, s, t, :]) : nothing
    end
    for (a, b) in zip(first, second)
        @assert size(a)[2:3] == size(b)[2:3]
        next = Nemo.zero_matrix(RR, size(a, 4), size(b, 4))
        for s in axes(a, 2), t in axes(a, 3)
            A, B = block(a, s, t), block(b, s, t)
            (isnothing(A) || isnothing(B)) && continue
            weight = isnothing(weights) ? F(1) : F(weights[t])/F(weights[s])
            next += RR(weight)*(transpose(A)*(G*B))
        end
        G = next
    end
    @assert Nemo.nrows(G) == Nemo.ncols(G) == 1
    return G[1, 1]
end

"Exact squared Hilbert--Schmidt norm of an MPO, using its ordinary double layer."
function hs_squared(cores; weights = nothing)
    G = Nemo.matrix(Nemo.QQ, ones(F, 1, 1))
    for core in cores
        l, dout, din, r = size(core)
        next = Nemo.zero_matrix(Nemo.QQ, r, r)
        for s in 1:dout, t in 1:din
            A = Nemo.matrix(Nemo.QQ, core[:, s, t, :])
            weight = isnothing(weights) ? F(1) : F(weights[t])/F(weights[s])
            next += weight*(transpose(A)*G*A)
        end
        G = next
    end
    @assert Nemo.nrows(G) == Nemo.ncols(G) == 1 && G[1, 1] >= 0
    return F(G[1, 1])
end

"Portable exact data for replay, with no dependency on the candidate generator."
function witness_encoding(input, certificate; precision = 192)
    data = Dict{String,Any}(
        "adjoint_first"=>hasproperty(certificate, :adjoint_first) && certificate.adjoint_first,
        "input"=>MN.encode_cores(input),
        "output"=>MN.encode_cores(certificate.cores),
        "bridges"=>[
            Dict("dimensions"=>collect(size(B)), "entries"=>MN.encode_fraction.(vec(B))) for
            B in certificate.bridges
        ],
        "residual_bounds"=>[string(r.residual_bound) for r in certificate.rows],
        "prefix_bounds"=>[string(r.local_prefix_bound) for r in certificate.rows],
        "operator_bounds"=>string.(certificate.tail.local_operator_bounds),
        "hs_bounds"=>string.(certificate.tail.local_hs_bounds),
        "tail_frobenius_capped"=>certificate.tail.frobenius_capped,
        "block_bounds"=>string.(certificate.tail.block_bounds),
        "error"=>string(certificate.error),
        "precision_bits"=>precision,
        "thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
        "error_method"=>string(
            hasproperty(certificate, :error_method) ? certificate.error_method : :triangle,
        ),
    )
    if get(data, "error_method", "triangle") == "projection"
        data["projection_bounds"] = [string(r.projection_bound) for r in certificate.rows]
    end
    if !isnothing(certificate.gauges)
        data["gauges"] = MN.encode_gauges(certificate.gauges)
    end
    if hasproperty(certificate, :weights) && !isnothing(certificate.weights)
        data["physical_weights"] = string.(certificate.weights)
    end
    if hasproperty(certificate, :operator_certificate) &&
       !isnothing(certificate.operator_certificate)
        independent = certificate.operator_certificate
        data["operator_certificate"] = Dict(
            "gauges"=>MN.encode_gauges(independent.gauges),
            "bounds"=>MN.encode_fraction.(independent.bounds),
        )
    end
    return data
end

"Read exact squaring/Gram data without invoking a candidate generator."
function decode_witness(data)
    input, output = MN.decode_cores(data["input"]), MN.decode_cores(data["output"])
    bridges = [
        reshape(MN.parse_fraction.(B["entries"]), Int.(B["dimensions"])...) for B in data["bridges"]
    ]
    rows = [
        (;
            residual_bound = MN.parse_fraction(data["residual_bounds"][j]),
            local_prefix_bound = MN.parse_fraction(data["prefix_bounds"][j]),
            projection_bound = haskey(data, "projection_bounds") ?
                               MN.parse_fraction(data["projection_bounds"][j]) : F(0),
        ) for j in eachindex(data["residual_bounds"])
    ]
    tail = (;
        local_operator_bounds = MN.parse_fraction.(data["operator_bounds"]),
        local_hs_bounds = MN.parse_fraction.(data["hs_bounds"]),
        frobenius_capped = get(data, "tail_frobenius_capped", false),
        block_bounds = MN.parse_fraction.(get(data, "block_bounds", String[])),
    )
    gauges = haskey(data, "gauges") ? MN.decode_gauges(data["gauges"]) : nothing
    weights =
        haskey(data, "physical_weights") ? MN.parse_fraction.(data["physical_weights"]) : nothing
    adjoint_first = get(data, "adjoint_first", false)
    operator_certificate = if haskey(data, "operator_certificate")
        independent = data["operator_certificate"]
        (;
            gauges = MN.decode_gauges(independent["gauges"]),
            bounds = MN.parse_fraction.(independent["bounds"]),
        )
    else
        nothing
    end
    certificate = (;
        cores = output,
        bridges,
        rows,
        tail,
        gauges,
        weights,
        adjoint_first,
        operator_certificate,
        error = MN.parse_fraction(data["error"]),
        error_method = Symbol(get(data, "error_method", "triangle")),
    )
    return input, certificate
end

"Decode a saved next candidate, or reconstruct its exact reflected-output encoding."
function flow_next_candidate(data; output_cores = nothing)
    if haskey(data, "next_candidate")
        return MN.decode_cores(data["next_candidate"])
    end
    @assert get(data, "next_candidate_encoding", nothing) == "reflected-output"
    @assert get(data, "reflected_output", false)
    kappa = MN.parse_fraction(data["normalization"])
    @assert kappa > 0
    C = isnothing(output_cores) ? MN.decode_cores(data["output"]) : output_cores
    next = MN.reflect_cores(C)
    next[1] = next[1]/kappa
    return next
end

"Operator norms certified backward from approximate Gram powers."
function backward_operator_bounds(normalizations, errors, terminal::F)
    @assert length(normalizations) == length(errors)
    @assert all(>(0), normalizations) && all(>=(0), errors) && terminal >= 0
    bounds = fill(F(0), length(errors)+1)
    bounds[end] = terminal
    for j in length(errors):-1:1
        bounds[j] = MN.relative_upper_sqrt(normalizations[j]*bounds[j+1]+errors[j])
    end
    return bounds
end

"Lipschitz bound for a Gram difference, using independently certified norms."
function gram_error_coefficient(error, candidate, reference = nothing)
    @assert error >= 0 && candidate >= 0
    isnothing(reference) && return 2candidate+error
    @assert reference >= 0
    # A'*A-B'*B = A'*(A-B)+(A'-B')*B. Both independent
    # norms can be much smaller than the available HS error bound.
    return min(2candidate+error, 2reference+error, candidate+reference)
end

"Forward Hilbert--Schmidt error propagation using already verified operator bounds."
function propagate_gram_errors(
    initial,
    normalizations,
    errors,
    bounds;
    round_bits = nothing,
    reference_bounds = nothing,
)
    @assert length(normalizations) == length(errors) && length(bounds) == length(errors)+1
    @assert isnothing(reference_bounds) || length(reference_bounds) == length(bounds)
    result = F[F(initial)]
    for j in eachindex(errors)
        e = last(result)
        reference = isnothing(reference_bounds) ? nothing : reference_bounds[j]
        coefficient = gram_error_coefficient(e, bounds[j], reference)
        next = (errors[j]+e*coefficient)/normalizations[j]
        push!(
            result,
            isnothing(round_bits) ? next : MN.relative_rational_upper(next; bits = round_bits),
        )
    end
    return result
end

end
