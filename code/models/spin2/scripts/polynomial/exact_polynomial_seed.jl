using Spin2Gap, JSON3, SHA
const PM = Spin2Gap.PolynomialMPO
const MN = PM.MN

function exact_seed_probe(output; lengths = (4, 16, 256), beta = 64, power = 28)
    result = Dict{String,Any}(
        "status"=>"exact-polynomial-MPO-compression",
        "complete"=>false,
        "uniform_gap_certified"=>false,
        "thermal_trace_certified"=>false,
        "runs"=>Any[],
    )
    for n in lengths
        started = time()
        seed = PM.PG.seed(n, beta, power)
        exact = PM.normalized_polynomial(seed)
        right = PM.exact_compress_right(exact)
        left = PM.exact_compress_right(MN.reflect_cores(right.cores))
        candidate = MN.reflect_cores(left.cores)
        proof = Dict(
            "length"=>n,
            "beta"=>string(seed.beta),
            "power"=>power,
            "remainder"=>string(seed.remainder),
            "exact_compression"=>true,
            "candidate"=>MN.encode_cores(candidate),
            "right_output"=>MN.encode_cores(right.cores),
            "right_transforms"=>PM.encode_transforms(right.transforms),
            "left_transforms"=>PM.encode_transforms(left.transforms),
            "error_upper"=>string(MN.F(0)),
            "precision_bits"=>192,
            "thermal_trace_certified"=>false,
            "uniform_gap_certified"=>false,
        )
        path = splitext(output)[1]*"-n$(n)-witness.json"
        open(path, "w") do io
            JSON3.write(io, proof)
        end
        dimensions = [size(c, 4) for c in candidate]
        push!(
            result["runs"],
            Dict(
                "length"=>n,
                "maxdim"=>maximum(dimensions),
                "dimensions"=>dimensions,
                "initial_error_upper"=>string(MN.F(0)),
                "witness_file"=>basename(path),
                "witness_sha256"=>bytes2hex(sha256(read(path))),
                "seconds"=>time()-started,
            ),
        )
        write_result(output, result)
        println(
            "Exact polynomial n=$n D=$(maximum(dimensions)) zero error, seconds=$(time()-started)",
        )
        flush(stdout)
    end
    result["complete"] = true
    write_result(output, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    exact_seed_probe(only(ARGS))
end
