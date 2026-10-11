# Exact fixed-point thermal contractions. The equations are ported from the
# Apache-2.0 upstream source1/source2; source2 uses separate integer coordinates
# instead of packing 192 lanes into a single giant integer.

const CENTERS1 = [
    124885379 65757;
    2925 2960742;
    12870614 286599;
    54069 1920414;
    4639283 510140;
    195455 1514846;
    2654823 675419;
    377629 1314544;
    1900547 787405
]
const CENTERS2 = [8945303 346734 939987; 3741281 567167 984108; 2327547 721090 1002613]

function poisson_coefficients(source)
    degree, mean, qbits = source == 1 ? (280, 98, 55) : (256, 84, 56)
    q = big(1) << qbits
    factorial_degree = factorial(big(degree))
    weights = [big(mean)^j * div(factorial_degree, factorial(big(j))) for j in 0:degree]
    total = sum(weights)
    coefficients = [Int64(fld(q * w, total)) for w in weights]
    while iszero(last(coefficients))
        pop!(coefficients)
    end
    @assert length(coefficients) == (source == 1 ? 192 : 173)
    @assert sum(big.(coefficients)) <= q
    return coefficients, q
end

function magnetization_basis(::Val{N}, magnetization) where {N}
    basis = NTuple{N,Int8}[]
    word = zeros(Int8, N)
    function visit(i, remaining)
        if i > N
            remaining == 0 && push!(basis, Tuple(word))
            return
        end
        for s in Int8(-1):Int8(1)
            abs(remaining - s) <= N - i || continue
            word[i] = s
            visit(i + 1, remaining - s)
        end
    end
    visit(1, magnetization)
    return basis
end

function representatives(basis::Vector{NTuple{N,Int8}}, source, magnetization) where {N}
    rows, multiplicities = Int[], Int[]
    flip = source == 1 && magnetization == 0
    for (row, word) in enumerate(basis)
        orbit = Set{NTuple{N,Int8}}()
        for reversed in (false, true), k in 0:(N-1)
            rotated = ntuple(i -> word[reversed ? N + 1 - mod1(i + k, N) : mod1(i + k, N)], N)
            push!(orbit, rotated)
            flip && push!(orbit, map(-, rotated))
        end
        if word == minimum(orbit)
            push!(rows, row)
            push!(multiplicities, length(orbit) * (magnetization == 0 ? 1 : 2))
        end
    end
    @assert sum(multiplicities) == length(basis) * (magnetization == 0 ? 1 : 2)
    return rows, multiplicities
end

# Each row lists (column, real coefficient, omega coefficient), omega^2=-1-omega.
function thermal_rows(basis::Vector{NTuple{N,Int8}}, source, twist) where {N}
    lookup = Dict(word => j for (j, word) in enumerate(basis))
    rows = Vector{Vector{Tuple{Int,Int64,Int64}}}(undef, length(basis))
    for (i, word) in enumerate(basis)
        bond_diagonal = sum(Int(word[j]) * word[mod1(j + 1, N)] for j in 1:N)
        diagonal = source == 1 ? 32 - 3N - 2bond_diagonal : 16 - div(3N, 2) - bond_diagonal
        row = [(i, Int64(diagonal), Int64(0))]
        for j in 1:N, delta in (-1, 1)
            k = mod1(j + 1, N)
            abs(word[j] - delta) <= 1 && abs(word[k] + delta) <= 1 || continue
            target = ntuple(l -> Int8(word[l] + (l == k ? delta : l == j ? -delta : 0)), N)
            a, b = source == 1 ? (-2, 0) : (-1, 0)
            if k == 1 && twist != 0
                a, b = twist == 1 ? (-a, 0) : delta == 1 ? (0, -1) : (1, 1)
            end
            push!(row, (lookup[target], a, b))
        end
        rows[i] = row
    end
    return rows
end

"Conservative integer bounds covering every partial sum, shift and addition."
function check_native_bounds(rows, source, n)
    _, q = poisson_coefficients(source)
    dimension_bound = isqrt(big(3)^n) + 1
    suffix_error = source == 1 ? 281 * (dimension_bound + 1) : big(125159)
    coordinate_bound = source == 1 ? q + suffix_error : 2 * (q + suffix_error)
    maximum_row = maximum(
        max(
            sum(abs(a) + abs(b) for (_, a, b) in row),
            sum(abs(b) + abs(a - b) for (_, a, b) in row),
        ) for row in rows
    )
    @assert maximum_row * coordinate_bound < typemax(Int64)
    @assert div(maximum_row * coordinate_bound, source == 1 ? 32 : 16) + q < typemax(Int64)
    return Int64(coordinate_bound)
end

function multiply_real_reference!(y, rows, x)
    lanes = size(x, 1)
    @inbounds for i in eachindex(rows)
        for lane in 1:lanes
            y[lane, i] = 0
        end
        for (j, a, _) in rows[i]
            @simd for lane in 1:lanes
                y[lane, i] += a * x[lane, j]
            end
        end
    end
end

function multiply_eisenstein_reference!(y, z, rows, u, v)
    lanes = size(u, 1)
    @inbounds for i in eachindex(rows)
        for lane in 1:lanes
            y[lane, i] = 0
            z[lane, i] = 0
        end
        for (j, a, b) in rows[i]
            @simd for lane in 1:lanes
                y[lane, i] += a * u[lane, j] - b * v[lane, j]
                z[lane, i] += b * u[lane, j] + (a - b) * v[lane, j]
            end
        end
    end
end

# At most 256 terms per Int128 chunk; this conservative bound fits even before
# the contraction bound is used. Checked adds also make overflow fail closed.
function squared_norm(u, v, lane, complex)
    total = big(0)
    chunk = Int128(0)
    for j in axes(u, 2)
        a = Int128(u[lane, j])
        b = complex ? Int128(v[lane, j]) : Int128(0)
        term = a * a - a * b + b * b
        chunk = Base.Checked.checked_add(chunk, term)
        if j % 256 == 0
            total += chunk
            chunk = 0
        end
    end
    return total + chunk
end

function thermal_sector(source, n, twist, magnetization; batch = 32, kernel = :optimized)
    kernel in (:reference, :optimized) || throw(ArgumentError("unknown integer kernel"))
    start = time()
    basis = magnetization_basis(Val(n), magnetization)
    reps, mult = representatives(basis, source, magnetization)
    rows = thermal_rows(basis, source, twist)
    bound = check_native_bounds(rows, source, n)
    coefficients, _ = poisson_coefficients(source)
    complex = source == 2 && twist == 2
    shift = source == 1 ? 5 : 4
    total = big(0)
    for first in 1:batch:length(reps)
        last = min(first + batch - 1, length(reps))
        rr = reps[first:last]
        u = zeros(Int64, length(rr), length(basis))
        y = similar(u)
        v = complex ? zero(u) : zeros(Int64, 0, 0)
        z = similar(v)
        for q in Iterators.reverse(coefficients)
            if complex
                if kernel == :optimized
                    multiply_eisenstein!(y, z, rows, u, v, shift)
                else
                    multiply_eisenstein_reference!(y, z, rows, u, v)
                    @. z = z >> shift
                    @. y = y >> shift
                end
                v, z = z, v
            elseif kernel == :optimized
                multiply_real!(y, rows, u, shift)
            else
                multiply_real_reference!(y, rows, u)
                @. y = y >> shift
            end
            for (lane, row) in enumerate(rr)
                y[lane, row] += q
            end
            u, y = y, u
        end
        @assert maximum(abs, u) <= bound
        complex && @assert maximum(abs, v) <= bound
        for (lane, multiplicity) in enumerate(mult[first:last])
            total += multiplicity * squared_norm(u, v, lane, complex)
        end
    end
    return (
        total = total,
        dimension = length(basis),
        representatives = length(reps),
        multiplicity = sum(mult),
        seconds = time() - start,
    )
end

function exp_interval(x, source)
    @assert 0 <= x < 15
    term = frac(1)
    lower = term
    for j in 1:100
        term *= x / j
        lower += term
    end
    return lower, lower + (source == 1 ? term : 2term * x / 101)
end

function thermal_certificate(source, n, twist, total)
    _, q = poisson_coefficients(source)
    a, beta = frac(700741, 500000), source == 1 ? frac(49, 4) : frac(21, 2)
    el, eu = exp_interval(beta * n * (frac(3, 2) - a), source)
    result =
        Dict{String,Any}("source" => source, "n" => n, "twist" => twist, "total" => string(total))
    if source == 1
        d = isqrt(big(3)^n) + 1
        err = d * (281 * (d + 1) + 1)
        r = isqrt(total)
        lower = el * frac(max(0, r - err)^2, q^2)
        upper = eu * frac((r + 1 + err)^2, q^2)
        center = frac(CENTERS1[n-3, twist+1], 10^6)
        epsilon = frac(30, 10^6)
        @assert center - epsilon < lower <= upper < center + epsilon
        @assert lower - (center - epsilon) > frac(13, 10^6)
        @assert center + epsilon - upper > frac(13, 10^6)
        expected = only(
            filter(
                x -> integer(x["n"]) == n && integer(x["tw"]) == twist,
                read_exact_json(joinpath(SOURCE, "evidence/source1/thermal-run/certificates.json")),
            ),
        )
        @assert total == integer(expected["W"])
        for (key, value) in
            ("lower" => lower, "upper" => upper, "center" => center, "allowed_error" => epsilon)
            @assert value == parse_rational(expected[key])
            result[key] = string(value)
        end
    else
        lower, upper = el * frac(total, q^2), eu * frac(total, q^2)
        center = frac(CENTERS2[div(n, 2)-2, twist+1], 10^6)
        @assert center - frac(2, 10^6) < lower <= upper < center + frac(2, 10^6)
        @assert upper < 16
        delta = 3^6 * 243 * (frac(257 * (2 * 243 + 1), q) + frac(1, big(2)^116))
        trace_error = (8 + delta) * delta
        @assert delta < frac(4, 10^7)
        @assert frac(2, 10^6) + trace_error < frac(1, 10^5)
        expected = read_exact_json(joinpath(SOURCE, "evidence/source2/runs/n$n-g$twist.json"))
        @assert total == integer(expected["total"])
        values = (
            "computed_trace_lower" => lower,
            "computed_trace_upper" => upper,
            "scaled_frobenius_error_bound" => delta,
            "trace_error_bound" => trace_error,
            "true_trace_lower" => lower - trace_error,
            "true_trace_upper" => upper + trace_error,
            "center" => center,
        )
        for (key, value) in values
            @assert value == parse_rational(expected[key])
            result[key] = string(value)
        end
    end
    result["matches_upstream_exactly"] = true
    return result
end

function run_thermal(
    source,
    output;
    lengths = source == 1 ? collect(4:12) : [6, 8, 10],
    batch = 32,
    kernel = :optimized,
)
    source in (1, 2) || throw(ArgumentError("source must be 1 or 2"))
    provenance = verify_source()
    records, sector_records = Any[], Any[]
    started = time()
    tasks = [(n, twist, m) for n in lengths for twist in (source == 1 ? (0:1) : (0:2)) for m in 0:n]
    # The upstream runner schedules all sectors together, largest chains first.
    # Greedy tasks also avoid leaving large sectors on one statically grouped thread.
    sort!(tasks; by = t -> (-t[1], t[3], t[2]))
    received = Dict{Tuple{Int,Int},Int}()
    totals = Dict{Tuple{Int,Int},BigInt}()
    multiplicities = Dict{Tuple{Int,Int},Int}()
    progress_lock = ReentrantLock()
    write_result(
        output,
        Dict(
            "status" => "running",
            "source_integrity" => provenance,
            "certificates" => records,
            "sectors" => sector_records,
        ),
    )
    Threads.@threads :greedy for task in tasks
        n, twist, m = task
        sector = thermal_sector(source, n, twist, m; batch, kernel)
        lock(progress_lock) do
            key = (n, twist)
            received[key] = get(received, key, 0) + 1
            totals[key] = get(totals, key, big(0)) + sector.total
            multiplicities[key] = get(multiplicities, key, 0) + sector.multiplicity
            push!(
                sector_records,
                Dict(
                    "n" => n,
                    "twist" => twist,
                    "magnetization" => m,
                    "total" => string(sector.total),
                    "dimension" => sector.dimension,
                    "representatives" => sector.representatives,
                    "multiplicity" => sector.multiplicity,
                    "wall_seconds" => sector.seconds,
                ),
            )
            if n >= 10 && m <= 3
                println(
                    "sector source=$source n=$n twist=$twist Sz=$m complete in $(round(sector.seconds; digits=2))s",
                )
            end
            if received[key] == n + 1
                @assert multiplicities[key] == 3^n
                push!(records, thermal_certificate(source, n, twist, totals[key]))
                println(
                    "PASS thermal source$source n=$n twist=$twist total=$(totals[key]) elapsed=$(round(time()-started; digits=1))s",
                )
                write_result(
                    output,
                    Dict(
                        "status" => "running",
                        "source_integrity" => provenance,
                        "certificates" => records,
                        "sectors" => sector_records,
                    ),
                )
            end
            flush(stdout)
        end
    end
    sort!(records; by = x -> (x["n"], x["twist"]))
    sort!(sector_records; by = x -> (x["n"], x["twist"], x["magnetization"]))
    return write_result(
        output,
        Dict(
            "status" => "reproduced",
            "method" => "fresh exact integer Julia contractions",
            "source_integrity" => provenance,
            "julia" => string(VERSION),
            "threads" => Threads.nthreads(),
            "batch" => batch,
            "kernel" => string(kernel),
            "schedule" => "largest-first global greedy sectors",
            "wall_seconds" => time() - started,
            "completed_utc" => string(now(UTC)),
            "certificates" => records,
            "sectors" => sector_records,
        ),
    )
end
