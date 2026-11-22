"Exact three-direction spacetime bootstrap, including intermediate square periods."
module CubeRecurrence

using CertifiedGaps.ExactArithmetic: RationalBig, frac, round_up, log_enclosure

export cube_step, cube_chain, direct_temperature_prefix

const GRID = big(10)^80
f(z) = z^2 / (2(1-z)^2)
q(e) = (1-e[1])^4 * (1-e[2])^2 * (1-e[3])

"One staggered x, y, thermal update; each output is rounded upward on a rational grid."
function cube_step(e; grid = GRID)
    x, y, t = RationalBig.(e)
    all(z -> 0 <= z < 1, (x, y, t)) || error("invalid input defects")
    grid > 0 || error("rounding grid must be positive")
    a = round_up(f(1 - q((x, y, t))), grid)
    a < 1 || error("x update leaves the admissible region")
    b = round_up(f(1 - q((y, t, a))), grid)
    b < 1 || error("y update leaves the admissible region")
    c = round_up(f(1 - q((t, a, b))), grid)
    c < 1 || error("thermal update leaves the admissible region")
    return (a, b, c)
end

"""
    direct_temperature_prefix(spatial, thermal; grid = 10^12)

Convert defects of Py(4,2,2β) and Pτ(4,4,β) into staggered cube defects at side four
and inverse temperature 2β. The face/cube identities and spatial/thermal squaring prove
this finite prefix; METHOD.md gives the indices and proof. No larger spectrum is needed.
"""
function direct_temperature_prefix(spatial, thermal; grid = big(10)^12)
    x, t = RationalBig(spatial), RationalBig(thermal)
    all(z -> 0 <= z < 1, (x, t)) || error("invalid direct-temperature input defects")
    grid > 0 || error("rounding grid must be positive")
    a = round_up(f(x), grid)
    a < 1 || error("direct spatial squaring leaves the admissible region")
    b = round_up(f(1 - (1-x)^2 * (1-a)), grid)
    b < 1 || error("direct face update leaves the admissible region")
    c = round_up(f(1 - (1-t)^4 * (1-a)^2 * (1-b)), grid)
    c < 1 || error("direct thermal update leaves the admissible region")
    return (a, b, c)
end

"""
    cube_chain(defects, beta; min_steps = 0, max_steps = 12, grid = GRID, min_even_side = 2)

Bound every even square torus of side at least `min_even_side`. At scale j, spatial
interpolation gives thermal purity at least Q_j on n₀2^j ≤ L ≤ n₀2^(j+1), where
n₀=`min_even_side`. For Q_j > 1/2 this yields gap ≥ log(Q_j/(1-Q_j))/(2^j beta).
Once max defect ≤ 1/1000, the scalar inequality e' ≤ 25e² gives the uniform tail
log(1/(25e))/(2^j beta). Return the minimum of all finite-window and tail bounds.
`min_steps` postpones the common scalar tail to retain unequal defects for longer.
The grid must resolve the smallest defects reached before that tail.
"""
function cube_chain(defects, beta; min_steps = 0, max_steps = 12, grid = GRID,
    min_even_side = 2)
    0 <= min_steps <= max_steps || error("invalid number of preliminary steps")
    grid > 0 || error("rounding grid must be positive")
    min_even_side isa Integer && min_even_side > 0 && iseven(min_even_side) ||
        error("the minimum side must be a positive even integer")
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
        if j >= min_steps && r <= frac(1, 1000)
            # This integer inequality proves the coefficient 25 on the entire tail interval.
            @assert big(49_000_000) < 50 * big(993)^2
            bound = first(log_enclosure(1 / (25r))) / b
            row["tail_gap_lower"] = string(bound)
            push!(bounds, bound)
            push!(rows, row)
            return Dict{String,Any}("rows" => rows,
                "gap_lower_bound" => string(minimum(bounds)),
                "rounding_grid" => string(grid), "min_even_side" => string(min_even_side))
        end
        Q = q(e)
        Q > 1//2 || error("intermediate-period purity is not above one half")
        bound = first(log_enclosure(Q / (1-Q))) / b
        row["cube_lower"] = string(Q)
        row["window_gap_lower"] = string(bound)
        push!(bounds, bound)
        push!(rows, row)
        e = cube_step(e; grid)
    end
    error("three-direction recurrence did not reach its quadratic tail")
end

end
