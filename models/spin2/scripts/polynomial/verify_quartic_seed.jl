using Spin2Gap, JSON3, SHA
const PM = Spin2Gap.PolynomialMPO
const MN = PM.MN
const MS = PM.MS

function verify_quartic_seed(path; residual_arithmetic = :interval)
    data = JSON3.read(read(path, String))
    @assert data.construction == "quartic-taylor-sos" && data.degree == 4
    read_factor = row -> begin
        witness = joinpath(dirname(path), String(row.witness_file))
        text = read(witness, String)
        @assert bytes2hex(sha256(codeunits(text))) == row.witness_sha256
        return JSON3.read(text)
    end
    progress = index -> begin
        println("Verified quartic factor $index/3")
        flush(stdout)
    end
    checked = PM.verify_quartic_seed_data(data, read_factor; residual_arithmetic, progress)
    println("Verified quartic seed at N=$(checked.seed.n): error<=$(Float64(checked.error))")
    return checked.error
end

if abspath(PROGRAM_FILE) == @__FILE__
    for path in ARGS
        verify_quartic_seed(path)
    end
end
