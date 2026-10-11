using Spin2Gap, ITensors, ITensorMPS, LinearAlgebra, JSON3, Dates, SHA
const ST = Spin2Gap.Trials
BLAS.set_num_threads(1)

function itensors_check(P, integer_check)
    sites = siteinds("S=2", 4)
    scale = Float64(maximum(maximum(abs, p) for p in P))
    A = [sqrt(ST.WEIGHTS[m])*Float64.(P[m])/scale for m in 1:5]
    amplitudes = [tr(A[6-a]*A[6-b]*A[6-c]*A[6-d]) for a in 1:5, b in 1:5, c in 1:5, d in 1:5]
    psi = MPS(ITensor(amplitudes, sites...), sites; cutoff = 0.0)
    normalize!(psi)
    H = Spin2Gap.hamiltonian(sites)
    energy = real(inner(psi', H, psi))
    exact = Float64(integer_check.E//integer_check.V)
    @assert abs(energy-exact)<1e-10
    os = OpSum()
    for j in 1:4
        os += 1.0, "S+", j
    end
    raising = MPO(os, sites)
    residual = sqrt(abs(real(inner(raising, psi, raising, psi))))
    @assert residual<1e-6
    return Dict(
        "length"=>4,
        "itensors_energy"=>energy,
        "exact_energy_display"=>exact,
        "absolute_energy_difference"=>abs(energy-exact),
        "total_raising_residual"=>residual,
        "certified"=>false,
    )
end

function main(output = joinpath(@__DIR__, "../results/trial-certificate.json"); witness = nothing)
    source = isnothing(witness) ? joinpath(@__DIR__, "../results/trial-polished.json") : witness
    input = JSON3.read(read(source, String), Dict{String,Any})
    input_hash = bytes2hex(sha256(read(source)))
    if isnothing(witness)
        row = only(filter(r->r["bond_dimension"]==24, input["runs"]))
        t = ST.templates(Int.(row["virtual_spins"]))
        parameters, P = ST.quantized_tensors(t, Float64.(row["x"]))
    else
        t = ST.templates(Int.(input["virtual_spins"]))
        parameters = parse.(BigInt, input["integer_parameters"])
        P = ST.polynomial_tensors(t, parameters)
    end
    Spin2Gap.write_result(output, Dict("status"=>"running", "uniform_gap_certified"=>false))
    started = time()
    blocks = ST.exact_blocks(t, P)
    println("Built exact trial blocks $(map(b->size(b.T,1),blocks))")
    flush(stdout)
    direct = ST.physical_integer_check(P)
    four = ST.exact_contract(blocks, 4)
    @assert four.V==direct.V && 4four.U==direct.E
    controls = itensors_check(P, direct)
    records = Any[]
    result = Dict(
        "status"=>"running",
        "certified_statement"=>"finite variational energy only",
        "uniform_gap_certified"=>false,
        "virtual_spins"=>t.spins,
        "bond_dimension"=>t.D,
        "integer_parameters"=>string.(parameters),
        "physical_weights"=>ST.WEIGHTS,
        "input_sha256"=>input_hash,
        "input_kind"=>isnothing(witness) ? "numerical-candidate-quantized-to-integers" :
                      "integer-witness-replay",
        "julia"=>string(VERSION),
        "Nemo"=>string(pkgversion(ST.Nemo)),
        "ITensors"=>string(pkgversion(ITensors)),
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "environment_sha256"=>Dict(
            name=>bytes2hex(sha256(read(joinpath(@__DIR__, "..", name)))) for
            name in ("Project.toml", "Manifest.toml")
        ),
        "itensors_control"=>controls,
        "certificates"=>records,
        "source_sha256"=>Dict(
            name=>bytes2hex(sha256(read(joinpath(@__DIR__, "../src", name)))) for
            name in ("Trials.jl", "trial_certificate.jl", "model.jl")
        ),
    )
    for n in (128, 256, 512)
        c = ST.exact_contract(blocks, n)
        if !isnothing(witness)
            expected = only(filter(row->row["length"]==n, input["certificates"]))
            @assert c.U==parse(BigInt, expected["U"]) && c.V==parse(BigInt, expected["V"])
            @assert c.records==expected["blocks"]
        end
        margin = -47595c.V-10000c.U
        @assert margin>0
        floor_energy = fld(big(10)^12*c.U, c.V)
        push!(
            records,
            Dict(
                "length"=>n,
                "U"=>string(c.U),
                "V"=>string(c.V),
                "energy_per_bond_floor_1e12"=>string(floor_energy),
                "shift_a"=>"9519//2000",
                "positive_margin_minus47595V_minus10000U"=>string(margin),
                "energy_per_bond_display"=>Float64(c.U//c.V),
                "blocks"=>c.records,
            ),
        )
        Spin2Gap.write_result(output, result)
        println("CERTIFIED N=$n: e=$(Float64(c.U//c.V)) < -4.7595 elapsed=$(time()-started)s")
        flush(stdout)
    end
    result["status"] = "finite-energy-certificates-passed"
    result["completed_utc"] = string(now(UTC))
    result["seconds"] = time()-started
    Spin2Gap.write_result(output, result)
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS)<=2 || error("usage: certify_trial.jl [OUTPUT.json] [INTEGER_WITNESS.json]")
    main(
        isempty(ARGS) ? joinpath(@__DIR__, "../results/trial-certificate.json") : ARGS[1];
        witness = length(ARGS)==2 ? ARGS[2] : nothing,
    )
end
