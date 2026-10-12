# Exact rational evaluation of the doubling chain with the bootstrap tail theorem of the
# periodic-gap source paper. Moved verbatim from spin1 `experiments/optimize_bounds.jl`
# (`gap_chain`) and `src/scalars.jl` (`purity_update`, `coefficient`); the output dictionary
# keys are unchanged because stored certificates use them.

export bootstrap_coefficient, rational_gap_chain, beta_window

"Bootstrap constant G(u) = (Σ_{j=1}^3 (1−u)^{−j})² / 2 of the uniform tail theorem."
bootstrap_coefficient(u) = sum((1 - u)^(-j) for j in 1:3)^2 / 2

# squaring_defect with the domain check spin1's `purity_update` made.
function _checked_squaring_defect(u)
    @assert 0 <= u < 1
    return squaring_defect(u)
end

"""
    rational_gap_chain(p, q, n, β; steps = 11, grid = 10^512, tail_C = nothing) -> Dict

Exact certificate of a uniform gap bound from rational seed defects p = 1 − Pₛ(n, β) and
q = 1 − P_τ(2n, β). It runs `steps` doublings of the exact map, rounding each new defect up
to the grid 1/`grid`, takes the per-scale bound log(A/(1−A))/β_j with A = (1−p)²(1−q), and
closes with the bootstrap tail: if u = max(p, q) at the last scale satisfies C·u < 1 with
C ≥ `bootstrap_coefficient(u)`, every later scale has gap ≥ log(1/(C u))/β. Logarithms are
enclosed by `log_enclosure`. Returns the record stored in the spin-1 certificates, with the
strict lower bound under `"gap_strict_lower"`.
"""
function rational_gap_chain(p::RationalBig, q::RationalBig, n::Int, beta::RationalBig;
    steps = 11, grid = big(10)^512, tail_C = nothing)
    initial_n = n
    initial_beta = beta
    intervals = Any[]
    # The running minimum starts at the first scale (spin1 started it at 1, which capped the
    # certified gap at 1; every spin-1 bound is far below that, so its records are unchanged).
    gap_lower = nothing
    for j in 0:(steps-1)
        @assert 0 < p < 1 && 0 < q < 1
        A = (1 - p)^2 * (1 - q)
        @assert frac(1, 2) < A < 1
        # The paper's interpolation gives T(L,beta) >= A for n <= L <= 2n.
        # Hence the ground state is unique and gamma_L >= log(A/(1-A))/beta.
        loglo, loghi = log_enclosure(A / (1 - A))
        gap_lower = isnothing(gap_lower) ? loglo / beta : min(gap_lower, loglo / beta)
        exact_p = _checked_squaring_defect(1 - (1 - p)^2 * (1 - q))
        next_p = round_up(exact_p, grid)
        exact_q = _checked_squaring_defect(1 - (1 - q)^2 * (1 - next_p))
        next_q = round_up(exact_q, grid)
        @assert exact_p <= next_p < exact_p + inv(frac(grid))
        @assert exact_q <= next_q < exact_q + inv(frac(grid))
        push!(
            intervals,
            Dict(
                "step" => j,
                "length_from" => n,
                "length_to" => 2n,
                "beta" => string(beta),
                "p" => string(p),
                "q" => string(q),
                "interpolated_purity_lower" => string(A),
                "gap_lower" => string(loglo / beta),
                "gap_upper_for_formula" => string(loghi / beta),
                "next_p" => string(next_p),
                "next_q" => string(next_q),
            ),
        )
        p, q = next_p, next_q
        n *= 2
        beta *= 2
    end
    u = max(p, q)
    G = bootstrap_coefficient(u)
    C = isnothing(tail_C) ? round_up(G, big(10)^24) : tail_C
    @assert 0 < u < frac(1, 6) && G <= C && C * u < 1 && C * (1 - 3u) > 3
    loglo, loghi = log_enclosure(inv(C * u))
    gap_lower = min(gap_lower, loglo / beta)
    # Strictly round below all interval bounds and the uniform tail theorem.
    claim = round_down(gap_lower, big(10)^9) - frac(1, big(10)^9)
    @assert 0 < claim < gap_lower
    return Dict(
        "even_length_threshold" => initial_n,
        "initial_beta" => string(initial_beta),
        "rounding_grid" => string(inv(frac(grid))),
        "gap_strict_lower" => string(claim),
        "gap_strict_lower_decimal" => Float64(claim),
        "finite_intervals" => intervals,
        "tail" => Dict(
            "length_threshold" => n,
            "beta" => string(beta),
            "p" => string(p),
            "q" => string(q),
            "u" => string(u),
            "C" => string(C),
            "G" => string(G),
            "r" => string(C * u),
            "gap_lower" => string(loglo / beta),
            "gap_upper_for_formula" => string(loghi / beta),
        ),
    )
end

"""
    beta_window(delta_energy, gap; purity = 0.99) -> Dict

Necessary single-scale window for a seed of the given `purity`, if `gap` and the spatial
defect `delta_energy` = E₀(2n) − 2E₀(n) were exact: β ≥ log(p/(1−p))/gap from the thermal
purity and β ≤ −2 log p / delta_energy from the spatial one. A screening heuristic, not a bound.
Moved verbatim from spin2 `scripts/screen.jl`.
"""
function beta_window(delta_energy, gap; purity = 0.99)
    0.5 < purity < 1 || throw(ArgumentError("purity must lie between 1/2 and 1"))
    gap > 0 || throw(ArgumentError("gap must be positive"))
    beta_min = log(purity / (1 - purity)) / gap
    beta_max = delta_energy > 0 ? -2log(purity) / delta_energy : nothing
    return Dict(
        "required_purity" => purity,
        "beta_min_from_thermal_purity" => beta_min,
        "beta_max_from_spatial_purity" => beta_max,
        "window_nonempty" => isnothing(beta_max) || beta_min <= beta_max,
        "delta_energy_over_gap" => delta_energy / gap,
        "necessary_ratio_ceiling" => -2log(purity) / log(purity / (1 - purity)),
        "certified" => false,
    )
end
