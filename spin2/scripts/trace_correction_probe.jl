# Reference experiment: retain the signed overlap in the final norm certificate.
# This certifies one compression, not a complete Gibbs approximation.
import ThermalMPOGap
using Spin2Gap, JSON3, SHA
import Nemo
const PM = Spin2Gap.PolynomialMPO
const MS = PM.MS
const MN = PM.MN
const MP = ThermalMPOGap.CertifiedMPO.Products

function trace_correction_probe(sources, output; reference_check = true)
    result = Dict{String,Any}(
        "status"=>"one-compression-trace-correction",
        "complete"=>false,
        "thermal_trace_certified"=>false,
        "uniform_gap_certified"=>false,
        "runs"=>Any[],
    )
    for source in sources
        started = time()
        B, proof = MS.decode_witness(JSON3.read(read(source, String)))
        @assert MS.verify_square_certificate(B, proof)
        weights = proof.weights
        first = proof.adjoint_first ? MN.adjoint_cores(B; weights) : B
        C = proof.cores
        cc = MS.hs_inner_ball(C, C; weights)
        ac = MP.product_inner_ball(C, first, B; weights)
        # The direct norm is used only as an independent reference for this small probe.
        RR = parent(cc)
        aa = if reference_check
            A=PM.product_mpo(first, B)
            MS.hs_inner_ball(A, A; weights)
        else
            nothing
        end
        corrected = 2PM.PG.relative_upper(ac)+PM.PG.relative_upper(-cc)+proof.error^2
        triangle = (PM.PG.relative_upper(sqrt(cc))+proof.error)^2
        @assert corrected >= 0
        if reference_check
            @assert aa-2ac+cc < RR(proof.error)^2
            @assert aa < RR(corrected) && aa < RR(triangle)
        end
        push!(
            result["runs"],
            Dict(
                "source_file"=>source,
                "source_sha256"=>bytes2hex(sha256(read(source))),
                "length"=>length(B),
                "max_output_dimension"=>maximum(size(c, 4) for c in C),
                "corrected_upper"=>string(corrected),
                "triangle_upper"=>string(triangle),
                "independent_full_reference_checked"=>reference_check,
                "reference_norm_squared_ball"=>isnothing(aa) ? nothing : string(aa),
                "overlap_ball"=>string(ac),
                "candidate_norm_squared_ball"=>string(cc),
                "certified_residual_squared_upper"=>string(proof.error^2),
                "triangle_over_reference_display"=>isnothing(aa) ? nothing :
                                                   Float64(RR(triangle)/aa),
                "corrected_over_reference_display"=>isnothing(aa) ? nothing :
                                                    Float64(RR(corrected)/aa),
                "triangle_over_candidate_display"=>Float64(RR(triangle)/cc),
                "corrected_over_candidate_display"=>Float64(RR(corrected)/cc),
                "improvement_over_candidate_norm_display"=>Float64(RR(triangle-corrected)/cc),
                "seconds"=>time()-started,
            ),
        )
        write_result(output, result)
        println(
            "$source: triangle/candidate=$(Float64(RR(triangle)/cc)), corrected/candidate=$(Float64(RR(corrected)/cc))",
        )
        flush(stdout)
    end
    result["complete"] = true
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) >= 2
    trace_correction_probe(ARGS[2:end], ARGS[1])
end
