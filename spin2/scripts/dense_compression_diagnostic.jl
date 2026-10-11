import ThermalMPOGap
using ThermalMPOGap.Results: write_result
using LinearAlgebra, JSON3, SHA
const MS = ThermalMPOGap.CertifiedMPO.Squaring
const MN = MS.MN
BLAS.set_num_threads(1)

"Floating-point four-site control, never a certificate input."
function dense_numeric(cores, weights)
    @assert length(cores) <= 4
    d = size(first(cores), 2)
    physical = sqrt.([Float64(weights[t]/weights[s]) for s in 1:d, t in 1:d])
    state = ones(1, 1, 1)
    for core in cores
        dim, _, left = size(state)
        right = size(core, 4)
        c = Float64.(core) .* reshape(physical, 1, d, d, 1)
        combined = reshape(state, dim*dim, left)*reshape(c, left, d*d*right)
        state = reshape(
            permutedims(reshape(combined, dim, dim, d, d, right), (3, 1, 4, 2, 5)),
            dim*d,
            dim*d,
            right,
        )
    end
    @assert size(state, 3) == 1
    return dropdims(state; dims = 3)
end

function dense_compression_diagnostic(sources, output)
    rows = Any[]
    for source in sources
        data = JSON3.read(read(source, String))
        witnesses = haskey(data, :runs) ? data.runs : [data]
        for (step, row) in enumerate(witnesses)
            path = joinpath(dirname(source), String(row.witness_file))
            text = read(path, String)
            @assert bytes2hex(sha256(codeunits(text))) == row.witness_sha256
            raw = JSON3.read(text)
            A, C, bound = if haskey(raw, :product)
                weights = MN.parse_fraction.(raw.product.weights)
                first = dense_numeric(MN.decode_cores(raw.product.first_input), weights)
                second = dense_numeric(MN.decode_cores(raw.product.second_input), weights)
                first*second,
                dense_numeric(MN.decode_cores(raw.product.output), weights),
                Float64(MN.parse_fraction(data.hs_error_upper))
            else
                input, certificate = MS.decode_witness(raw)
                B = dense_numeric(input, certificate.weights)
                (certificate.adjoint_first ? transpose(B)*B : B*B),
                dense_numeric(certificate.cores, certificate.weights),
                Float64(certificate.error)
            end
            E = A-C
            result = Dict(
                "source"=>source,
                "step"=>step,
                "hs_error_estimate"=>norm(E),
                "operator_error_estimate"=>opnorm(E),
                "saved_hs_error_bound"=>bound,
                "relative_hs_error_estimate"=>norm(E)/norm(A),
                "floating_scale_indicator"=>eps(Float64)*(norm(A)+norm(C)),
            )
            push!(rows, result)
            println(result)
            flush(stdout)
        end
    end
    write_result(
        output,
        Dict(
            "certified"=>false,
            "status"=>"floating-point-small-system-diagnostic",
            "uniform_gap_certified"=>false,
            "rows"=>rows,
        ),
    )
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) >= 2
    dense_compression_diagnostic(ARGS[1:(end-1)], last(ARGS))
end
