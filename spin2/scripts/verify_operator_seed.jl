import ThermalMPOGap
using Spin2Gap, JSON3, SHA
const PM = Spin2Gap.PolynomialMPO
const OE = ThermalMPOGap.CertifiedMPO.OperatorError
const MP = ThermalMPOGap.CertifiedMPO.Products
const MN = PM.MN

function checked_exact_seed(path)
    text = read(path, String)
    data = JSON3.read(text)
    @assert get(data, :exact_compression, false)
    seed = PM.scalar_seed(data)
    @assert seed.remainder == MN.parse_fraction(data.remainder)
    @assert MN.parse_fraction(data.error_upper) == 0
    candidate = MN.decode_cores(data.candidate)
    @assert PM.verify_exact_seed(PM.normalized_polynomial(seed), candidate, data)
    return (; seed, candidate, hash = bytes2hex(sha256(codeunits(text))))
end

"Replay an exact polynomial source and its compressed operator-norm certificate."
function verify_operator_seed(path; residual_arithmetic = :interval)
    data = JSON3.read(read(path, String))
    source = checked_exact_seed(String(data.source_file))
    @assert source.hash == data.source_sha256
    witness_path = joinpath(dirname(path), String(data.witness_file))
    text = read(witness_path, String)
    @assert bytes2hex(sha256(codeunits(text))) == data.witness_sha256
    raw = JSON3.read(text)
    A, identity, product = MP.decode(raw.product)
    @assert A == source.candidate && identity == PM.identity_mpo(source.seed.n)
    @assert product.weights == PM.WEIGHTS
    @assert MP.verify_product(A, identity, product; residual_arithmetic)
    operator = OE.decode(raw.operator)
    @assert OE.verify_operator_error(
        A,
        product,
        operator;
        second_input = identity,
        residual_arithmetic,
    )
    @assert MN.parse_fraction(data.operator_error_upper) >= operator.error
    @assert MN.parse_fraction(data.hs_error_upper) >= product.error
    @assert MN.decode_cores(data.candidate) == product.cores
    @assert Int(data.length) == source.seed.n
    println(
        "Verified compressed polynomial seed at N=$(source.seed.n): operator error <= $(Float64(operator.error))",
    )
    return (;
        candidate = product.cores,
        operator_error = operator.error,
        hs_error = product.error,
        seed = source.seed,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    foreach(verify_operator_seed, ARGS)
end
