"""
    ThermalMPOGap.CertifiedMPO.PolynomialBounds

Exact Taylor thermal majorants Q(K) ≥ exp(-δK) with a scalar inflation bound, Arb-accepted
rational upper endpoints, and the positive sum-of-squares form of even Taylor polynomials.
Moved verbatim from spin2 `Spin2PolynomialBounds`, except that the spectral radius of
K = H + nI, hard-coded there as 5n (the spin-2 bond spectrum), is now the keyword `radius`,
whose default `5n` reproduces spin2 exactly.
"""
module PolynomialBounds

import Nemo
using ..CertifiedMPO: F

"Accept a rational upper endpoint by an outward-rounded ball comparison."
function rational_upper(x; bits = 128)
    denominator = big(2)^bits
    center = round(BigInt, BigFloat(x)*BigFloat(denominator))
    radius = big(2)
    RR = parent(x)
    while !(x < RR((center+radius)//denominator))
        radius *= 2
    end
    return (center+radius)//denominator
end

"A rational ball upper endpoint with relative dyadic resolution."
function relative_upper(x; bits = 128)
    value = BigFloat(x)
    exponent2 = iszero(value) ? 0 : exponent(abs(value))
    scale = exponent2 >= 0 ? F(big(2)^exponent2) : F(1, big(2)^(-exponent2))
    return scale*rational_upper(x/parent(x)(scale); bits)
end

"""
Exact polynomial thermal majorant, requiring no Trotter formula.

For K=H+nI, the spin-2 bond spectrum gives -5n I <= K <= 5n I.
Other models pass their own radius M with -M I <= K <= M I (default M = 5n).
Let delta=beta/(2m), m=2^power, and choose rational
r >= exp(delta*M)*(delta*M)^3/6. Then
Q=I-delta*K+delta^2*K^2/2+r*I >= exp(-delta*K)>0.
Both are functions of K, so eigenvalue comparison gives
Z_H(beta) <= exp(beta*n)*Tr(Q^(2m)).

The scalar inflation bound compares this EXACT polynomial trace with Z.
It does not certify an MPO approximation or contraction of that trace.
"""
function seed(n, beta, power; precision = 192, radius = 5n)
    @assert n isa Integer && power isa Integer
    @assert iseven(n) && n >= 4 && beta > 0 && power >= 0
    @assert beta isa Union{Integer,Rational}
    RR = Nemo.ArbField(precision)
    m = big(2)^power
    delta = F(beta)/(2m)
    M = radius
    x = RR(delta*M)
    remainder = rational_upper(exp(x)*x^3/6)
    # |Q-e^{-delta K}| <= 2r, and e^{-delta K} >= e^{-delta M}I.
    # Thus Q^(2m) <= (1+2r e^{delta M})^(2m) e^{-beta K}.
    # 1+y <= e^y supplies the recorded simpler scalar bound.
    log_inflation = 4m*RR(remainder)*exp(x)
    inflation_upper = rational_upper(exp(log_inflation))
    return (;
        n,
        beta = F(beta),
        m,
        power,
        delta,
        spectral_radius = M,
        remainder,
        identity = F(1)+remainder,
        linear = -delta,
        quadratic = delta^2/2,
        inflation_upper,
        log_inflation,
        precision,
    )
end

polynomial(c, x) = c.identity+c.linear*x+c.quadratic*x^2

"""
Higher-order Taylor majorant with the same spectral comparison as `seed`.

For degree d and x=delta*radius (default radius 5n), choose r >= exp(x)*x^(d+1)/(d+1)!.
Then Q(K)=sum((-delta*K)^k/k!, k=0:d)+r*I >= exp(-delta*K)>0.
The scalar log-inflation is at most 4m*r*exp(x). This certifies only
the exact polynomial; a compressed MPO still needs its own error witness.
The existing quadratic seed and its saved witnesses are unchanged.
"""
function taylor_seed(n, beta, power, degree; precision = 192, radius = 5n)
    @assert n isa Integer && power isa Integer && degree isa Integer
    @assert iseven(n) && n >= 4 && beta > 0 && power >= 0 && degree >= 1
    @assert beta isa Union{Integer,Rational}
    RR = Nemo.ArbField(precision)
    m = big(2)^power
    delta = F(beta)/(2m)
    x = RR(delta*radius)
    remainder = relative_upper(exp(x)*x^(degree+1)/factorial(big(degree+1)))
    coefficients = F[F(1)+remainder]
    coefficient = F(1)
    for k in 1:degree
        coefficient *= -delta/k
        push!(coefficients, coefficient)
    end
    log_inflation_upper = relative_upper(4m*RR(remainder)*exp(x))
    operator_upper = sum(abs(coefficients[k+1])*big(radius)^k for k in 0:degree)
    return (;
        n,
        beta = F(beta),
        power,
        degree,
        m,
        delta,
        spectral_radius = radius,
        remainder,
        coefficients,
        operator_upper,
        log_inflation_upper,
        precision,
    )
end

function taylor_polynomial(c, x)
    value = last(c.coefficients)
    for coefficient in Iterators.reverse(c.coefficients[1:(end-1)])
        value = value*x+coefficient
    end
    return value
end

"""
Exact positive sum of squares for the even Taylor polynomial P_(2k)(x).

P_(2k)(x) = integral_0^infinity exp(-t)*(t+x)^(2k) dt / (2k)!.
The moment Gram matrix of (t+x)^k is positive definite. Rational LDL'
therefore expresses P_(2k) as positive rational weights times squares of
degree-k polynomials. Coefficient equality is checked exactly below.
"""
function even_taylor_sos(degree)
    @assert degree isa Integer && degree >= 2 && iseven(degree)
    k = div(degree, 2)
    G = [
        F(binomial(big(k), i)*binomial(big(k), j)*factorial(big(2k-i-j)), factorial(big(2k)))
        for i in 0:k, j in 0:k
    ]
    L, weights = zeros(F, k+1, k+1), zeros(F, k+1)
    for j in 1:(k+1)
        L[j, j] = 1
        weights[j] = G[j, j]-sum((L[j, h]^2*weights[h] for h in 1:(j-1)); init = F(0))
        @assert weights[j] > 0
        for i in (j+1):(k+1)
            L[i, j] =
                (G[i, j]-sum((L[i, h]*L[j, h]*weights[h] for h in 1:(j-1)); init = F(0)))/weights[j]
        end
    end
    polynomials = [L[:, j] for j in 1:(k+1)]
    coefficients = zeros(F, degree+1)
    for (weight, polynomial) in zip(weights, polynomials), i in 0:k, j in 0:k
        coefficients[i+j+1] += weight*polynomial[i+1]*polynomial[j+1]
    end
    @assert coefficients == [F(1, factorial(big(j))) for j in 0:degree]
    return (; weights, polynomials)
end

end
