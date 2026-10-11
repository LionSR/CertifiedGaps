import ThermalMPOGap
using Spin2Gap, LinearAlgebra, JSON3, SHA, Dates
import Nemo
const TB = Spin2Gap.ThermalBounds
const PB = ThermalMPOGap.CertifiedMPO.PositiveBounds
const F = PB.F
BLAS.set_num_threads(1)

function amplitude_interval(beta, m; shift = 238061//50000)
    @assert iseven(m) && m >= 2 && beta > 0
    n=4
    RR = Nemo.ArbField(192)
    bond = TB.rational_bond()
    halfgate = TB.gate_matrix(bond, RR(beta)/(2m), RR)
    gate = TB.gate_matrix(bond, RR(beta)/m, RR)
    M = zeros(F, 625, 625)
    lower = zeros(F, 625, 625)
    exact_trace = RR(0)
    powers = [5^(j-1) for j in 1:n]
    prefactor = exp(-RR(beta)*RR(shift)*n/2)
    for mag in -8:8
        states = [code for code in 0:624 if sum(mod(div(code, p), 5)-2 for p in powers) == mag]
        a = TB.layer_matrix(n, mag, halfgate, 1, RR)
        b = TB.layer_matrix(n, mag, gate, 2, RR)
        Q = a*b*a
        amplitude = prefactor*Q^div(m, 2)
        exact_trace += prefactor^2*Nemo.tr(Q^m)
        for (row, r) in enumerate(states), (col, c) in enumerate(states)
            # Reshape the operator into its Schmidt matrix across sites 12|34.
            left = mod(r, 25)+25mod(c, 25)+1
            right = div(r, 25)+25div(c, 25)+1
            lo, upper = PB.rational_bracket(amplitude[row, col])
            M[left, right] = max(upper, zero(F))
            lower[left, right] = max(lo, zero(F))
        end
    end
    return lower, M, exact_trace
end

function amplitude_upper(beta, m; kwargs...)
    _, upper, trace = amplitude_interval(beta, m; kwargs...)
    return upper, trace
end

function local_charge(code)
    ket, bra = mod(code, 25), div(code, 25)
    return mod(ket, 5)+div(ket, 5)-mod(bra, 5)-div(bra, 5)
end

function local_metric(code)
    ket, bra = mod(code, 25), div(code, 25)
    weight(a) = binomial(4, mod(a, 5))*binomial(4, div(a, 5))
    return F(weight(bra), weight(ket))
end

function main(
    output = joinpath(@__DIR__, "../results/positive-compression-probe.json");
    rankcaps = (4, 8, 16),
    iterations = 1500,
    underweight = 256.0,
    cover_sweeps = 0,
    correction_rank = 1,
    betas = (1, 4),
)
    result = Dict{String,Any}(
        "status"=>"four-site-positive-compression-probe",
        "complete"=>false,
        "uniform_gap_certified"=>false,
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "nmf_iterations"=>iterations,
        "underestimation_penalty"=>underweight,
        "cover_sweeps"=>cover_sweeps,
        "correction_rank"=>correction_rank,
        "runs"=>Any[],
        "source_sha256"=>Dict(
            name=>bytes2hex(sha256(read(joinpath(@__DIR__, "../src", name)))) for
            name in ("../../julia/src/CertifiedMPO/PositiveBounds.jl", "ThermalBounds.jl")
        ),
    )
    metric = [local_metric(i) for i in 0:624]
    groups = Dict(q=>findall(i->local_charge(i-1)==q, 1:625) for q in -8:8)
    for beta in betas
        started = time()
        M, trace_bound = amplitude_upper(beta, 64)
        initial_upper = sum(M[i, j]^2*metric[i]*metric[j] for i in 1:625, j in 1:625)
        RR = parent(trace_bound)
        @assert trace_bound < initial_upper < trace_bound*(1+RR(1)/big(10)^18)
        for rankcap in rankcaps
            bound = zeros(F, 625, 625)
            totalrank = 0
            for q in -8:8
                rows, cols = groups[q], groups[-q]
                c = PB.upper_factorization(
                    M[rows, cols],
                    min(rankcap, length(rows), length(cols));
                    iterations,
                    row_metric = metric[rows],
                    col_metric = metric[cols],
                    underweight,
                    cover_sweeps,
                    correction_rank,
                )
                bound[rows, cols] = c.bound
                totalrank += c.rank
            end
            @assert all(bound .>= M)
            # Frobenius norm after returning to the orthonormal magnetic basis.
            upperZ = sum(bound[i, j]^2*metric[i]*metric[j] for i in 1:625, j in 1:625)
            @assert trace_bound < upperZ
            ratio = Float64(upperZ)/Float64(trace_bound)
            row = Dict(
                "length"=>4,
                "beta"=>beta,
                "steps"=>64,
                "rankcap_per_charge"=>rankcap,
                "total_bond_dimension"=>totalrank,
                "entrywise_majorant_verified_exactly"=>true,
                "thermal_trace_upper_rational"=>string(upperZ),
                "trotter_trace_ball"=>string(trace_bound),
                "compression_factor_display"=>ratio,
                "seconds"=>time()-started,
            )
            push!(result["runs"], row)
            write_result(output, result)
            println(
                "positive compression beta=$beta rankcap=$rankcap D=$totalrank upper/exact-Trotter=$ratio elapsed=$(time()-started)s",
            )
            flush(stdout)
        end
    end
    result["complete"] = true
    result["completed_utc"] = string(now(UTC))
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(
        isempty(ARGS) ? joinpath(@__DIR__, "../results/positive-compression-probe.json") :
        only(ARGS),
    )
end
