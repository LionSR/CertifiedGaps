using Spin2Gap, JSON3, SHA, Dates
import Nemo
const SC = Spin2Gap.Certificate
const ST = SC.ST

function enclosure(energy; denominator = big(10)^12)
    center = round(BigInt, Float64(energy)*Float64(denominator))
    radius = big(2)
    while !(
        energy > Nemo.QQ(center-radius, denominator) &&
        energy < Nemo.QQ(center+radius, denominator)
    )
        radius *= 2
        radius <= denominator || error("energy enclosure too broad")
    end
    return Dict(
        "lower_numerator"=>string(center-radius),
        "upper_numerator"=>string(center+radius),
        "denominator"=>string(denominator),
    )
end

function main(
    output,
    source;
    dimension = 89,
    replay = false,
    lengths = (128, 256, 512),
    scale = 10^8,
    verify_images = true,
    split_rotation = true,
    flint_threads = get(ENV, "NEMO_THREADED", "0") == "1" ? min(4, Sys.CPU_THREADS) : 1,
)
    started = time()
    flint_threads > 1 && Nemo.flint_set_num_threads(flint_threads)
    input = JSON3.read(read(source, String), Dict{String,Any})
    if replay
        t = ST.templates(Int.(input["virtual_spins"]))
        parameters = parse.(BigInt, input["integer_parameters"])
    else
        candidate = only(filter(r -> r["bond_dimension"] == dimension, input["runs"]))
        t = ST.templates(Int.(candidate["virtual_spins"]))
        parameters, _ = ST.quantized_tensors(t, Float64.(candidate["x"]); scale)
    end
    P = ST.polynomial_tensors(t, parameters)
    records = Any[]
    result = Dict{String,Any}(
        "status"=>"running",
        "complete"=>false,
        "certified_statement"=>"finite variational energy only",
        "uniform_gap_certified"=>false,
        "virtual_spins"=>t.spins,
        "bond_dimension"=>t.D,
        "integer_parameters"=>string.(parameters),
        "source_sha256"=>bytes2hex(sha256(read(source))),
        "input_kind"=>replay ? "integer-witness-replay" : "quantized-numerical-candidate",
        "implementation_sha256"=>Dict(
            name=>bytes2hex(sha256(read(joinpath(@__DIR__, "../src", name)))) for
            name in ("Certificate.jl", "Trials.jl", "trial_certificate.jl")
        ),
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "environment_sha256"=>Dict(
            name=>bytes2hex(sha256(read(joinpath(@__DIR__, "..", name)))) for
            name in ("Project.toml", "Manifest.toml")
        ),
        "julia"=>string(VERSION),
        "Nemo"=>string(pkgversion(Nemo)),
        "flint_threads"=>flint_threads,
        "certificates"=>records,
    )
    function save()
        result["seconds"] = time()-started
        result["updated_utc"] = string(now(UTC))
        write_result(output, result)
    end
    save()
    reduced = SC.reduced_blocks(t, P; verify_images, progress = true)
    result["multiplicity_dimensions"] = length.(getproperty.(reduced, :pairs))
    result["highest_weight_images_checked_exactly"] = verify_images
    result["tensor_covariance_checked_exactly"] = true
    result["complete_irreducible_bases_checked_exactly"] = true
    result["reduction_justification"] = "Exact covariance and complete irreducible bases; Schur lemma determines each image from one pivot coefficient."
    result["reduction_seconds"] = time()-started
    save()
    if split_rotation
        reduced = SC.split_transpose_rotation(reduced; progress = true)
    end
    result["transpose_rotation_commutation_checked_exactly"] = split_rotation
    result["arithmetic_block_dimensions"] = length.(getproperty.(reduced, :pairs))
    save()
    precision = 192
    contractions = SC.ball_contract_many(reduced, lengths; precision, scale, progress = true)
    for n in lengths
        c = contractions[n]
        interval = enclosure(c.energy)
        # A strict, exactly checked rational threshold is reported with each result.
        a =
            -Nemo.QQ(
                parse(BigInt, interval["upper_numerator"]),
                parse(BigInt, interval["denominator"]),
            )
        @assert c.energy < -a
        push!(
            records,
            Dict(
                "length"=>n,
                "energy_per_bond_interval"=>interval,
                "energy_per_bond_display"=>Float64(c.energy),
                "ball_display"=>string(c.energy),
                "precision_bits"=>precision,
                "shift_a"=>string(a),
                "strict_variational_inequality_verified"=>true,
            ),
        )
        save()
        println("CERTIFIED finite trial D=$(t.D) N=$n e=$(c.energy) elapsed=$(time()-started)s")
        flush(stdout)
    end
    result["status"] = "finite-energy-certificates-passed"
    result["complete"] = true
    save()
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        main(
            joinpath(@__DIR__, "../results/trial-reduced-certificate.json"),
            joinpath(@__DIR__, "../results/trial-analytic.json"),
        )
    else
        length(ARGS) in (3, 4) ||
            error("usage: certify_reduced_trial.jl OUTPUT INPUT DIMENSION [replay]")
        main(
            ARGS[1],
            ARGS[2];
            dimension = parse(Int, ARGS[3]),
            replay = length(ARGS)==4 && ARGS[4]=="replay",
        )
    end
end
