"""
    Spin2Gap.PolynomialMPO

The spin-2-specific polynomial thermal seeds: the centered Hamiltonian K = H + N I in the
Marshall basis, the degree-two and quartic (Taylor SOS) majorant MPOs, and their exact replay.
The d-generic MPO algebra (direct sums, products, exact compression and conditioning, initial
errors, sparse expansion) lives in `ThermalMPOGap.CertifiedMPO.MPOAlgebra`; its names are
re-bound here so that `PolynomialMPO.linear_combination` etc. keep working.
"""
module PolynomialMPO

using LinearAlgebra, SparseArrays
import Nemo
import ThermalMPOGap.CertifiedMPO: Squaring, Norms, PolynomialBounds, MPOAlgebra
using ThermalMPOGap.CertifiedMPO.MPOAlgebra:
    WEIGHTS,
    identity_mpo,
    linear_combination,
    product_core,
    product_mpo,
    exact_compress_right,
    verify_exact_compression,
    encode_transforms,
    decode_transforms,
    condition_right_exact,
    initial_error,
    sparse_operator
const MS = Squaring
const MN = Norms
const PG = PolynomialBounds
const MA = MPOAlgebra
const F = MN.F

function local_spins()
    plus, minus = zeros(F, 5, 5), zeros(F, 5, 5)
    for k in 0:3
        plus[k+2, k+1] = 4-k
        minus[k+1, k+2] = k+1
    end
    return (; identity = Matrix{F}(I, 5, 5), z = Matrix(Diagonal(F.(-2:2))), plus, minus)
end

"K=H+N I for the even periodic spin-2 chain, exactly in the Marshall basis."
function centered_hamiltonian(n)
    @assert iseven(n) && n >= 4
    s = local_spins()
    W = zeros(F, 5, 5, 5, 5)
    W[1, :, :, 1] = W[5, :, :, 5] = s.identity
    W[1, :, :, 5] = s.identity # Sum one identity per site, giving N I.
    W[1, :, :, 2], W[2, :, :, 5] = s.z, s.z
    W[1, :, :, 3], W[3, :, :, 5] = s.plus, -s.minus/2
    W[1, :, :, 4], W[4, :, :, 5] = s.minus, -s.plus/2
    open_chain = vcat([W[1:1, :, :, :]], fill(W, n-2), [W[:, :, :, 5:5]])
    wraps = Vector{Array{F,4}}[]
    for (a, b) in ((s.z, s.z), (-s.plus/2, s.minus), (-s.minus/2, s.plus))
        push!(wraps, vcat([reshape(a, 1, 5, 5, 1)], identity_mpo(n-2), [reshape(b, 1, 5, 5, 1)]))
    end
    return linear_combination(vcat([open_chain], wraps), ones(F, 4))
end

"The exact degree-two thermal majorant from its certified scalar coefficients."
function polynomial_mpo(seed)
    K = centered_hamiltonian(seed.n)
    return linear_combination(
        [identity_mpo(seed.n), K, product_mpo(K, K)],
        [seed.identity, seed.linear, seed.quadratic],
    )
end

function normalized_polynomial(seed)
    Q = polynomial_mpo(seed)
    Q[1] = Q[1]/big(5)^div(seed.n, 2)
    return Q
end

"Norm upper bound for the exact normalized seed, from the two scalar endpoints."
function initial_operator_upper(seed)
    if hasproperty(seed, :coefficients)
        return seed.operator_upper/big(5)^div(seed.n, 2)
    end
    radius = F(5seed.n)
    return (seed.identity+abs(seed.linear)*radius+seed.quadratic*radius^2)/big(5)^div(seed.n, 2)
end

"Reconstruct the exact scalar polynomial specified by a saved seed."
function scalar_seed(data)
    n, beta, power = Int(data["length"]), MN.parse_fraction(data["beta"]), Int(data["power"])
    degree = Int(get(data, "degree", 2))
    if degree == 2
        return PG.seed(n, beta, power)
    end
    @assert degree == 4 && get(data, "construction", nothing) == "quartic-taylor-sos"
    return PG.taylor_seed(n, beta, power, degree)
end

"Norm bounds for exact polynomial powers under the recorded positive rescalings."
function reference_operator_bounds(seed, kappas; precision = 192)
    RR = Nemo.ArbField(precision)
    bounds = F[initial_operator_upper(seed)]
    for kappa in kappas
        @assert kappa > 0
        next = PG.relative_upper(RR(last(bounds))^2/RR(kappa))
        @assert next*kappa >= last(bounds)^2
        push!(bounds, next)
    end
    return bounds
end

"The three exact quadratic factors of the fourth-order Taylor polynomial."
function quartic_sos_factors(seed)
    @assert seed.degree == 4
    sos = PG.even_taylor_sos(4)
    factors = [
        (; n = seed.n, identity = p[1], linear = -seed.delta*p[2], quadratic = seed.delta^2*p[3]) for p in sos.polynomials
    ]
    return (; sos.weights, factors)
end

"Certify one normalized quadratic-factor Gram approximation."
function verify_quartic_factor(
    seed,
    index,
    input,
    certificate;
    precision = 192,
    residual_arithmetic = :interval,
    progress = nothing,
    exact_input = nothing,
)
    @assert certificate.adjoint_first && certificate.weights == WEIGHTS
    factor = quartic_sos_factors(seed).factors[index]
    exact = normalized_polynomial(factor)
    initial_epsilon = if isnothing(exact_input)
        initial_error(exact, input; precision).error
    else
        @assert verify_exact_seed(exact, input, exact_input)
        F(0)
    end
    checked = MS.verify_square_certificate(
        input,
        certificate;
        details = true,
        residual_arithmetic,
        progress,
    )
    a = initial_operator_upper(factor)
    b = min(checked.operator_norm_upper, checked.hs_norm_upper)
    coefficient = MS.gram_error_coefficient(initial_epsilon, b, a)
    error = MN.relative_rational_upper(certificate.error+coefficient*initial_epsilon)
    return (; error, initial_error = initial_epsilon, local_square_error = certificate.error)
end

"Assemble Q/5^(N/2) from normalized quadratic-factor Gram approximations."
function quartic_sos_assembly(seed, outputs)
    @assert seed.degree == 4 && length(outputs) == 3
    sigma = F(1, big(5)^div(seed.n, 2))
    components = [copy(C) for C in outputs]
    for C in components
        C[1] = C[1]/sigma
    end
    identity = identity_mpo(seed.n)
    identity[1] = sigma*identity[1]
    return linear_combination(
        vcat(components, [identity]),
        vcat(quartic_sos_factors(seed).weights, [seed.remainder]),
    )
end

"Propagate the three verified Gram errors and the final compression residual."
function quartic_sos_seed_error(seed, outputs, errors, candidate; precision = 192)
    @assert length(errors) == 3 && all(>=(0), errors)
    combined = quartic_sos_assembly(seed, outputs)
    compression = initial_error(combined, candidate; precision)
    sigma = F(1, big(5)^div(seed.n, 2))
    factor_error = sum(quartic_sos_factors(seed).weights .* errors)/sigma
    error = MN.relative_rational_upper(factor_error+compression.error)
    return (; error, factor_error, compression_error = compression.error)
end

"Replay all quartic seed inputs; the caller supplies hash-checked witness data."
function verify_quartic_seed_data(
    data,
    read_factor;
    residual_arithmetic = :interval,
    progress = nothing,
)
    @assert data["construction"] == "quartic-taylor-sos" && data["degree"] == 4
    seed = scalar_seed(data)
    @assert seed.remainder == MN.parse_fraction(data["remainder"])
    @assert length(data["factors"]) == 3
    precision = Int(data["precision_bits"])
    outputs, errors = Vector{Array{F,4}}[], F[]
    for (index, row) in enumerate(data["factors"])
        @assert row["index"] == index
        raw = read_factor(row)
        input, certificate = MS.decode_witness(raw)
        checked = verify_quartic_factor(
            seed,
            index,
            input,
            certificate;
            precision,
            residual_arithmetic,
            exact_input = get(raw, "exact_input", nothing),
        )
        @assert MN.parse_fraction(row["gram_error_upper"]) >= checked.error
        push!(outputs, certificate.cores)
        push!(errors, checked.error)
        isnothing(progress) || progress(index)
    end
    candidate = MN.decode_cores(data["candidate"])
    checked = quartic_sos_seed_error(seed, outputs, errors, candidate; precision)
    claimed = MN.parse_fraction(data["error_upper"])
    @assert claimed >= checked.error
    return (; seed, candidate, error = claimed)
end

"Replay a two-sweep exact polynomial compression from portable local identities."
function verify_exact_seed(exact, candidate, data)
    if get(data, "exact_conditioning", false)
        unconditioned = MN.decode_cores(data["unconditioned_candidate"])
        @assert verify_exact_compression(
            unconditioned,
            candidate,
            decode_transforms(data["conditioning_transforms"]),
        )
        candidate = unconditioned
    end
    right = MN.decode_cores(data["right_output"])
    @assert verify_exact_compression(exact, right, decode_transforms(data["right_transforms"]))
    @assert verify_exact_compression(
        MN.reflect_cores(right),
        MN.reflect_cores(candidate),
        decode_transforms(data["left_transforms"]),
    )
    return true
end

end
