"Exact evaluations of the graph-level Ising high-field bound proved in METHOD.md."
module HighField

using CertifiedGaps, JSON3
using CertifiedGaps.ExactArithmetic: RationalBig, parse_rational
using CertifiedGaps.Results

export field_bound, certify_high_field, check_high_field, FIELD_POINTS

const FIELD_POINTS = RationalBig[41//10, 17//4, 9//2, 5, 6, 8, 18, 20, 24]
const FORMAT = "CertifiedGaps Ising path bound v1"
const GRID = big(10)^12

"Enclose 2√(h(h-s)) by exact rational arithmetic, for minimum field h > incident sum s."
function field_bound(h, s; grid = GRID)
    h, s = RationalBig(h), RationalBig(s)
    h > s >= 0 || error("the graph theorem requires h > s >= 0")
    grid isa Integer && grid > 0 || error("the square-root grid must be a positive integer")
    grid = big(grid)
    radicand = 4h*(h-s)
    scaled = radicand*grid^2
    p = isqrt(div(numerator(scaled), denominator(scaled)))
    lower = p//grid
    upper = p^2 == scaled ? lower : (p+1)//grid
    @assert 0 <= lower && lower^2 <= radicand <= upper^2
    return Dict("minimum_field" => string(h), "incident_sum_bound" => string(s),
        "gap_squared" => string(radicand), "gap_lower" => string(lower),
        "gap_upper_for_formula" => string(upper), "rounding_grid" => string(grid))
end

record_path(root) = joinpath(root, "results", "high_field_bounds.json")

function payload()
    return Dict("format" => FORMAT,
        "method" => "ferromagnetic Ising path expansion and FKG covariance bound",
        "scope" => "finite graphs with nonnegative pair couplings, g_i >= h > s >= sum_j J_ij",
        "square_lattice" => "J = 1; every finite open or periodic rectangle; g > 4",
        "formalization" => "mathematical proof in METHOD.md; not a Lean model theorem",
        "points" => [field_bound(g, 4) for g in FIELD_POINTS])
end

"Write the exact arithmetic record; its physical justification is the theorem in METHOD.md."
function certify_high_field(root)
    record = payload()
    record["metadata"] = parse_exact_json(JSON3.write(model_metadata(root,
        Dict("method" => "Ising path bound"); dirs = ("src", "scripts"),
        packages = (CertifiedGaps,), certified = true)))
    write_result(record_path(root), record)
    return record
end

"Check every exact number and the source hashes of the analytic-bound record."
function check_high_field(root)
    record = read_result(record_path(root))
    verify_model_metadata(record["metadata"], root; dirs = ("src", "scripts"))
    expected = payload()
    Set(keys(record)) == union(Set(keys(expected)), Set(["metadata"])) ||
        error("unexpected analytic-bound record fields")
    all(record[k] == v for (k, v) in expected) || error("analytic-bound replay mismatch")
    return record
end

end
