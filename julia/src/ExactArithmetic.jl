"""
    CertifiedGaps.ExactArithmetic

Exact rational arithmetic with outward rounding. Every routine returns rationals, and every
enclosure is valid without floating point, so its output can certify an inequality.

Moved verbatim from spin1 (`optimize_bounds.jl`, `temperature_lift.jl`) and spin2
(`Spin2MPONorms.jl`); only the rounding scale became an explicit argument.
"""
module ExactArithmetic

export RationalBig, frac, round_up, round_down, log_enclosure, first_true, integer_root,
    power_interval, power_upper, encode_fraction, parse_fraction, parse_rational

"Exact rationals with arbitrary-precision numerator and denominator."
const RationalBig = Rational{BigInt}

"The rational a/b with `BigInt` numerator and denominator."
frac(a::Integer, b::Integer = 1) = BigInt(a) // BigInt(b)

"Smallest point of the grid (1/scale)ℤ that is ≥ x."
round_up(x, scale) = frac(cld(numerator(x) * scale, denominator(x)), scale)

"Largest point of the grid (1/scale)ℤ that is ≤ x."
round_down(x, scale) = frac(fld(numerator(x) * scale, denominator(x)), scale)

"""
    log_enclosure(x) -> (lower, upper)

Outward rational enclosure of log x for rational x ≥ 1, by range reduction to [1, 2) and the
atanh series with its remainder. Both ends lie on the grid 10⁻²⁴.
"""
function log_enclosure(x::RationalBig)
    @assert x >= 1
    k = ndigits(numerator(x); base = 2) - ndigits(denominator(x); base = 2)
    k = max(k, 0)
    while frac(2)^k > x
        k -= 1
    end
    while frac(2)^(k + 1) <= x
        k += 1
    end
    y = x / frac(2)^k
    @assert 1 <= y < 2
    function series(z)
        @assert 0 <= z <= frac(1, 3)
        terms = 40
        partial = 2sum(z^(2j + 1) / (2j + 1) for j in 0:(terms-1))
        remainder = 2z^(2terms + 1) / ((2terms + 1) * (1 - z^2))
        return partial, partial + remainder
    end
    # Rounding the mantissa first keeps the rational powers small.
    grid = big(10)^24
    yl, yu = round_down(y, grid), round_up(y, grid)
    lower, _ = series((yl - 1) / (yl + 1))
    _, upper = series((yu - 1) / (yu + 1))
    l2, u2 = series(frac(1, 3))
    return round_down(k * l2 + lower, grid), round_up(k * u2 + upper, grid)
end

"""
    first_true(predicate, lo, hi; digits = 12) -> (a, b)

Bisection on the decimal grid 10^-digits for an increasing predicate: returns adjacent grid
points with `predicate(a) == false` and `predicate(b) == true`.
"""
function first_true(predicate, lo, hi; digits = 12)
    scale = big(10)^digits
    a, b = fld(numerator(lo) * scale, denominator(lo)), cld(numerator(hi) * scale, denominator(hi))
    @assert !predicate(frac(a, scale)) && predicate(frac(b, scale))
    while b - a > 1
        middle = fld(a + b, 2)
        if predicate(frac(middle, scale))
            b = middle
        else
            a = middle
        end
    end
    @assert !predicate(frac(a, scale)) && predicate(frac(b, scale))
    return frac(a, scale), frac(b, scale)
end

"Floor of the positive integer bth root, without floating-point arithmetic."
function integer_root(n::BigInt, b::Int)
    @assert n >= 0 && b >= 1
    b == 1 && return n
    b == 2 && return isqrt(n)
    n == 0 && return big(0)
    lo = big(0)
    hi = big(1) << cld(ndigits(n; base = 2), b)
    while hi - lo > 1
        mid = (lo + hi) >> 1
        if mid^b <= n
            lo = mid
        else
            hi = mid
        end
    end
    @assert lo^b <= n < (lo + 1)^b
    return lo
end

"Outward enclosure of x^r for nonnegative rational x and r, on the grid 10^-digits."
function power_interval(x::RationalBig, r::RationalBig; digits = 24)
    @assert x >= 0 && r >= 0
    r == 0 && return (frac(1), frac(1))
    x == 0 && return (frac(0), frac(0))
    a, b = Int(numerator(r)), Int(denominator(r))
    b == 1 && return (x^a, x^a)
    scale = big(10)^digits
    target = x^a
    n = fld(numerator(target) * scale^b, denominator(target))
    root = integer_root(n, b)
    lo = frac(root, scale)
    hi = lo^b == target ? lo : frac(root + 1, scale)
    @assert lo^b <= target <= hi^b
    return lo, hi
end

"Upper end of `power_interval`."
power_upper(x, r) = last(power_interval(x, r))

"""
    encode_fraction(x) -> String

Exact token for a rational: `"num//den"`, or `"num//2^k"` when the denominator is a power of
two and the exponent form is shorter.
"""
function encode_fraction(x::Rational)
    den = denominator(x)
    if ispow2(den)
        exponent = ndigits(den; base = 2) - 1
        token = "2^" * string(exponent)
        if length(token) < ndigits(den)
            return string(numerator(x)) * "//" * token
        end
    end
    return string(x)
end

"Inverse of `encode_fraction`; requires the `num//den` or `num//2^k` form."
function parse_fraction(s)
    parts = split(s, "//")
    @assert length(parts) == 2
    den = if startswith(parts[2], "2^")
        exponent = parse(Int, parts[2][3:end])
        @assert exponent >= 0
        big(1) << exponent
    else
        parse(BigInt, parts[2])
    end
    @assert !iszero(den)
    return parse(BigInt, parts[1]) // den
end

"""
    parse_rational(x)

Exact rational from an integer, an integer string, `"a/b"`, `"a//b"`, or a dictionary with
`numerator` and `denominator`. Refuses floating-point input.
"""
parse_rational(x::Integer) = frac(x)
parse_rational(x::AbstractDict) =
    frac(_integer(x["numerator"]), _integer(x["denominator"]))
function parse_rational(x::AbstractString)
    parts = split(x, r"/+")
    return length(parts) == 1 ? frac(_integer(parts[1])) :
           frac(_integer(parts[1]), _integer(parts[2]))
end
parse_rational(x) = error("Refusing inexact rational input: $(typeof(x))")

_integer(x::AbstractString) = parse(BigInt, x)
_integer(x::Integer) = BigInt(x)
_integer(x) = error("Refusing inexact integer input: $(typeof(x))")

end
