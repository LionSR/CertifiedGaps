using Spin2Gap, JSON3, SHA
const PM = Spin2Gap.PolynomialMPO
const MN = PM.MN

function verify_seed(path; residual_arithmetic = :interval)
    data = JSON3.read(read(path, String))
    seed = PM.scalar_seed(data)
    @assert seed.remainder == MN.parse_fraction(data.remainder)
    if get(data, :construction, nothing) == "quartic-taylor-sos"
        read_factor = row -> begin
            witness = joinpath(dirname(path), String(row.witness_file))
            text = read(witness, String)
            @assert bytes2hex(sha256(codeunits(text))) == row.witness_sha256
            return JSON3.read(text)
        end
        progress = index -> begin
            println("Verified source quartic factor $index/3")
            flush(stdout)
        end
        checked = PM.verify_quartic_seed_data(data, read_factor; residual_arithmetic, progress)
        println("Verified quartic polynomial seed error <= $(Float64(checked.error)): $path")
        return checked.error
    end
    exact = PM.normalized_polynomial(seed)
    candidate = MN.decode_cores(data.candidate)
    claimed = MN.parse_fraction(data.error_upper)
    if get(data, :exact_compression, false)
        @assert PM.verify_exact_seed(exact, candidate, data)
        @assert claimed == 0
    else
        checked = PM.initial_error(exact, candidate; precision = Int(data.precision_bits))
        @assert claimed >= checked.error
    end
    println("Verified normalized polynomial seed error <= $(Float64(claimed)): $path")
    return claimed
end

if abspath(PROGRAM_FILE) == @__FILE__
    for path in ARGS
        verify_seed(path)
    end
end
