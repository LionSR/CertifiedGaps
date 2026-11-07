"Exact three-direction spacetime bootstrap, including intermediate square periods."
module CubeRecurrence

using CertifiedGaps.ExactArithmetic: RationalBig, frac, round_up, log_enclosure

export cube_step, cube_chain

const GRID = big(10)^80
f(z) = z^2 / (2(1-z)^2)
q(e) = (1-e[1])^4 * (1-e[2])^2 * (1-e[3])

"One staggered x, y, thermal update; each output is rounded upward on a rational grid."
function cube_step(e)
    x, y, t = RationalBig.(e)
    all(z -> 0 <= z < 1, (x, y, t)) || error("invalid input defects")
    a = round_up(f(1 - (1-x)^4 * (1-y)^2 * (1-t)), GRID)
    a < 1 || error("x update leaves the admissible region")
    b = round_up(f(1 - (1-y)^4 * (1-t)^2 * (1-a)), GRID)
    b < 1 || error("y update leaves the admissible region")
    c = round_up(f(1 - (1-t)^4 * (1-a)^2 * (1-b)), GRID)
    c < 1 || error("thermal update leaves the admissible region")
    return (a, b, c)
end

"""
    cube_chain(defects, beta; max_steps = 12)

Bound all even square tori of side at least two. At scale j the interpolation in the
2+1 dimensional note gives thermal purity at least Q_j on 2^(j+1) ≤ L ≤ 2^(j+2).
For Q_j > 1/2 this yields gap ≥ log(Q_j/(1-Q_j))/(2^j beta). Once max defect ≤ 1/1000,
the scalar inequality e' ≤ 25e² gives the uniform tail log(1/(25e))/(2^j beta).
The returned rational lower bound is the minimum of the finite windows and the tail.
"""
function cube_chain(defects, beta; max_steps = 12)
    beta = RationalBig(beta)
    beta > 0 || error("inverse temperature must be positive")
    e = Tuple(RationalBig.(defects))
    all(z -> 0 < z < 1, e) || error("defects must lie strictly between zero and one")
    rows = Dict{String,Any}[]
    bounds = RationalBig[]
    for j in 0:max_steps
        b = beta * big(2)^j
        r = maximum(e)
        row = Dict{String,Any}("scale" => string(j), "beta" => string(b),
            "defects" => string.(collect(e)))
        if r <= frac(1, 1000)
            # This integer inequality proves the coefficient 25 on the entire tail interval.
            @assert big(49_000_000) < 50 * big(993)^2
            bound = first(log_enclosure(1 / (25r))) / b
            row["tail_gap_lower"] = string(bound)
            push!(bounds, bound)
            push!(rows, row)
            return Dict{String,Any}("rows" => rows,
                "gap_lower_bound" => string(minimum(bounds)),
                "rounding_grid" => string(GRID), "min_even_side" => "2")
        end
        Q = q(e)
        Q > 1//2 || error("intermediate-period purity is not above one half")
        bound = first(log_enclosure(Q / (1-Q))) / b
        row["cube_lower"] = string(Q)
        row["window_gap_lower"] = string(bound)
        push!(bounds, bound)
        push!(rows, row)
        e = cube_step(e)
    end
    error("three-direction recurrence did not reach its quadratic tail")
end

end
