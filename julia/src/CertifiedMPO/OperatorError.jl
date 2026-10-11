"""
    CertifiedGaps.CertifiedMPO.OperatorError

Telescoping operator-norm error of a compressed square or product, with independent vector or
matrix gauges for each of the three residual families, streamed residuals, replay and witness
encoding. Moved verbatim from spin2 `Spin2MPOOperatorError`; `MS` is `CertifiedMPO.Squaring`.
"""
module OperatorError

using ..CertifiedMPO: nemo_array
import ..Squaring
const MS = Squaring
const MN = MS.MN
const F = MN.F
const Nemo = MS.Nemo
using LinearAlgebra

unit_gauges(cores) = vcat([F[1]], [ones(F, size(c, 4)) for c in cores])

first_input(input, square) =
    hasproperty(square, :adjoint_first) && square.adjoint_first ?
    MN.adjoint_cores(input; weights = square.weights) : input

input_gauges(square) = hasproperty(square, :gauges) ? square.gauges : square.first_gauges

"Validate vector or matrix gauges and compute exact inverse identities."
function gauge_matrices(cores, gauges)
    @assert length(gauges) == length(cores)+1
    dimensions = vcat(size(first(cores), 1), [size(c, 4) for c in cores])
    matrices, inverses = Matrix{F}[], Matrix{F}[]
    for (g, d) in zip(gauges, dimensions)
        if g isa AbstractVector
            @assert length(g) == d && all(>(0), g)
            matrix = Matrix(Diagonal(F.(g)))
        else
            @assert g isa AbstractMatrix && size(g) == (d, d)
            matrix = F.(g)
        end
        exact = Nemo.matrix(Nemo.QQ, matrix)
        inverse = inv(exact)
        @assert exact*inverse == Nemo.identity_matrix(Nemo.QQ, d)
        push!(matrices, matrix)
        push!(inverses, nemo_array(inverse))
    end
    @assert first(matrices) == last(matrices) == ones(F, 1, 1)
    return (; matrices, inverses)
end

"Change all three residual auxiliary bases without forming a Kronecker matrix."
function transformed_residual(input, square, j, old_core, prefix, first_right, second_right; field)
    bridge = j == 1 ? ones(F, 1, 1, 1) : square.bridges[j-1]
    Z = MS.square_bridge(bridge, old_core, input[j]; field)
    chi, d, di, r, rb = size(Z)
    @assert d == di
    C, next = square.cores[j], square.bridges[j]
    @assert size(C)[1:3] == (chi, d, d) && size(next) == (size(C, 4), r, rb)
    E =
        reshape(Z, chi*d*d, r*rb)-MS.matrix_product(
            reshape(C, chi*d*d, size(C, 4)),
            reshape(next, size(C, 4), r*rb),
            field,
        )
    E = reshape(E, chi, d, d, r, rb)
    left = reshape(MS.matrix_product(prefix, reshape(E, chi, d*d*r*rb), field), chi, d, d, r, rb)
    # In Julia's column order the first tail index precedes the second one.
    first_axis = reshape(permutedims(left, (1, 2, 3, 5, 4)), chi*d*d*rb, r)
    middle = MS.matrix_product(first_axis, first_right, field)
    middle = permutedims(reshape(middle, chi, d, d, rb, r), (1, 2, 3, 5, 4))
    right = MS.matrix_product(reshape(middle, chi*d*d*r, rb), second_right, field)
    return reshape(right, chi, d, d, r*rb)
end

"Visit local residuals one at a time; the callback may request an exact fallback."
function foreach_residual(
    callback,
    input,
    square,
    first_gauges,
    second_gauges,
    prefix_gauges;
    arithmetic = :rational,
    precision = 192,
    second_input = input,
)
    @assert arithmetic in (:rational, :interval)
    n = length(input)
    @assert n == length(square.cores) == length(second_input) && length(square.bridges) == n-1
    weights = square.weights
    first = first_input(input, square)
    old_gauges = isnothing(input_gauges(square)) ? unit_gauges(input) : input_gauges(square)
    old = gauge_matrices(first, old_gauges)
    first_data = gauge_matrices(first, first_gauges)
    second_data = gauge_matrices(second_input, second_gauges)
    prefix_data = gauge_matrices(square.cores, prefix_gauges)
    original_first = MN.gauge_transform(
        first,
        old_gauges;
        interval_precision = arithmetic == :interval ? precision : nothing,
    )
    field = arithmetic == :interval ? Nemo.ArbField(precision) : Nemo.QQ
    for j in 1:(n-1)
        first_right = MS.rational_product(old.matrices[j+1], first_data.inverses[j+1])
        second_right, prefix = second_data.inverses[j+1], prefix_data.matrices[j]
        e = transformed_residual(
            second_input,
            square,
            j,
            original_first[j],
            prefix,
            first_right,
            second_right;
            field,
        )
        exact =
            () ->
                arithmetic == :rational ? e :
                transformed_residual(
                    second_input,
                    square,
                    j,
                    MN.exact_gauge_core(first[j], old_gauges[j], old_gauges[j+1]),
                    prefix,
                    first_right,
                    second_right;
                    field = Nemo.QQ,
                )
        callback(j, e, field, exact)
    end
    # No truncation at the final core: its equality closes the telescoping sum.
    bridge = n == 1 ? ones(F, 1, 1, 1) : square.bridges[end]
    terminal = MS.square_bridge(bridge, original_first[end], second_input[end])
    C = last(square.cores)
    @assert size(C, 4) == 1 && size(terminal)[4:5] == (1, 1)
    @assert size(C)[1:3] == size(terminal)[1:3] && vec(C) == vec(terminal)
    return nothing
end

"Collect exact residuals for small controls; certificate generation streams them."
function residuals(input, square, first_gauges, second_gauges, prefix_gauges; second_input = input)
    result = Array{F,4}[]
    foreach_residual(
        input,
        square,
        first_gauges,
        second_gauges,
        prefix_gauges;
        second_input,
    ) do _, e, _, _
        push!(result, e)
    end
    return result
end

"The exact telescoping sum, with linear-time prefix and tail products."
function telescoping_error(bounds, residual_bounds)
    n = length(bounds[1])
    @assert length(bounds) == 3 && all(length(b) == n for b in bounds)
    @assert length(residual_bounds) == n-1
    tails = ones(F, n+1)
    for j in n:-1:1
        tails[j] = tails[j+1]*bounds[1][j]*bounds[2][j]
    end
    prefix, error = F(1), F(0)
    for j in eachindex(residual_bounds)
        error += prefix*residual_bounds[j]*tails[j+1]
        prefix *= bounds[3][j]
    end
    return error
end

"Operator-error bound; all proposed local norms are checked using ball arithmetic."
function operator_error_certificate(
    input,
    square;
    second_input = input,
    first_gauges = unit_gauges(input),
    second_gauges = unit_gauges(second_input),
    prefix_gauges = unit_gauges(square.cores),
    precision = 192,
    arithmetic = :rational,
    progress = nothing,
)
    weights = square.weights
    first = first_input(input, square)
    bounds = Vector{F}[]
    for (cores, gauges) in
        ((first, first_gauges), (second_input, second_gauges), (square.cores, prefix_gauges))
        maps = MN.gauge_transform(
            cores,
            gauges;
            interval_precision = arithmetic == :interval ? precision : nothing,
        )
        push!(bounds, [MN.core_norm_upper(c, :operator; weights, precision) for c in maps])
    end
    residual_bounds = F[]
    foreach_residual(
        input,
        square,
        first_gauges,
        second_gauges,
        prefix_gauges;
        arithmetic,
        precision,
        second_input,
    ) do j, e, field, _
        metric = MN.hs_row_metric(size(e, 1), size(e, 2), weights)
        push!(
            residual_bounds,
            MS.frobenius_upper(reshape(e, length(metric), size(e, 4)), field; metric),
        )
        isnothing(progress) ||
            progress((; stage = :operator_residual, site = j, total = length(input)))
    end
    error = telescoping_error(bounds, residual_bounds)
    return (; first_gauges, second_gauges, prefix_gauges, bounds, residual_bounds, error)
end

"Replay local inequalities and the telescoping operator-error bound without an optimizer."
function verify_operator_error(
    input,
    square,
    proof;
    precision = 192,
    residual_arithmetic = :rational,
    progress = nothing,
    second_input = input,
)
    @assert residual_arithmetic in (:rational, :interval)
    weights = square.weights
    first = first_input(input, square)
    @assert length(proof.bounds) == 3
    families = (
        (first, proof.first_gauges),
        (second_input, proof.second_gauges),
        (square.cores, proof.prefix_gauges),
    )
    for (k, (cores, gauges)) in enumerate(families)
        maps = MN.gauge_transform(
            cores,
            gauges;
            interval_precision = residual_arithmetic == :interval ? precision : nothing,
        )
        @assert length(proof.bounds[k]) == length(input)
        for (c, alpha) in zip(maps, proof.bounds[k])
            @assert MN.verify_core_bound(c, :operator, alpha; weights, precision)
        end
    end
    @assert length(proof.residual_bounds) == length(input)-1
    foreach_residual(
        input,
        square,
        proof.first_gauges,
        proof.second_gauges,
        proof.prefix_gauges;
        arithmetic = residual_arithmetic,
        precision,
        second_input,
    ) do j, e, field, exact
        alpha = proof.residual_bounds[j]
        metric = MN.hs_row_metric(size(e, 1), size(e, 2), weights)
        @assert alpha >= 0
        squared = sum(
            abs2.(reshape(e, length(metric), size(e, 4))) .*
            (field === Nemo.QQ ? metric : field.(metric)),
        )
        if field === Nemo.QQ
            @assert alpha^2 >= squared
        elseif !(field(alpha)^2 > squared)
            rational = exact()
            @assert alpha^2 >=
                    sum(abs2.(reshape(rational, length(metric), size(rational, 4))) .* metric)
        end
        isnothing(progress) ||
            progress((; stage = :operator_residual, site = j, total = length(input)))
    end
    error = telescoping_error(proof.bounds, proof.residual_bounds)
    @assert error == proof.error
    return true
end

function encode(proof)
    return Dict(
        "first_gauges"=>MN.encode_gauges(proof.first_gauges),
        "second_gauges"=>MN.encode_gauges(proof.second_gauges),
        "prefix_gauges"=>MN.encode_gauges(proof.prefix_gauges),
        "bounds"=>[string.(b) for b in proof.bounds],
        "residual_bounds"=>string.(proof.residual_bounds),
        "error"=>string(proof.error),
    )
end

function decode(data)
    vectors(key) = [MN.parse_fraction.(x) for x in data[key]]
    return (;
        first_gauges = MN.decode_gauges(data["first_gauges"]),
        second_gauges = MN.decode_gauges(data["second_gauges"]),
        prefix_gauges = MN.decode_gauges(data["prefix_gauges"]),
        bounds = vectors("bounds"),
        residual_bounds = MN.parse_fraction.(data["residual_bounds"]),
        error = MN.parse_fraction(data["error"]),
    )
end

end
