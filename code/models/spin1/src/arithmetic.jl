# Independent construction from COLUMN actions, versus the production ROW rules.
eadd(a, b) = (a[1] + b[1], a[2] + b[2])
emul(a, b) = (a[1] * b[1] - a[2] * b[2], a[1] * b[2] + a[2] * b[1] - a[2] * b[2])
econj(a) = (a[1] - a[2], -a[2])
phase(g, k) =
    g == 0 ? (1, 0) : g == 1 ? (iseven(k) ? 1 : -1, 0) : ((1, 0), (0, 1), (-1, -1))[mod(k, 3)+1]

function physical_hamiltonian(::Val{N}, twist) where {N}
    words = [Tuple(w) for w in Iterators.product(ntuple(_ -> -1:1, N)...)]
    h = Dict{Tuple{NTuple{N,Int},NTuple{N,Int}},Tuple{Int,Int}}()
    for y in words
        h[(y, y)] = (sum(y[j] * y[mod1(j + 1, N)] for j in 1:N), 0)
        for j in 1:N, step in (-1, 1)
            k = mod1(j + 1, N)
            x = ntuple(l -> y[l] + (l == j ? step : l == k ? -step : 0), N)
            maximum(abs, x) <= 1 || continue
            value = k == 1 ? phase(twist, y[1] - x[1]) : (1, 0)
            h[(x, y)] = eadd(get(h, (x, y), (0, 0)), value)
        end
    end
    return h
end

function matrix_and_symmetry_checks()
    records = Any[]
    for source in (1, 2),
        n in (source == 1 ? (4, 5, 6) : (4, 6)),
        twist in (source == 1 ? (0:1) : (0:2))

        h = physical_hamiltonian(Val(n), twist)
        represented_entries = 0
        for magnetization in (-n):n
            basis = magnetization_basis(Val(n), magnetization)
            rows = thermal_rows(basis, source, twist)
            for (i, row) in enumerate(rows), (j, a, b) in row
                x, y = Int.(basis[i]), Int.(basis[j])
                coefficient = source == 1 ? -2 : -1
                offset = source == 1 ? 32 - 3n : 16 - div(3n, 2)
                expected = eadd(emul((coefficient, 0), h[(x, y)]), x == y ? (offset, 0) : (0, 0))
                @assert (a, b) == expected
                represented_entries += 1
            end
        end
        @assert represented_entries == length(h)
        for ((x, y), value) in h
            @assert value == econj(h[(y, x)]) && sum(x) == sum(y)
            @assert h[(reverse(x), reverse(y))] == econj(value)
            @assert h[(map(-, x), map(-, y))] == econj(value)
            sx, sy = (x[end], x[1:(end-1)]...), (y[end], y[1:(end-1)]...)
            @assert h[(sx, sy)] == emul(emul(phase(twist, -x[end]), value), phase(twist, y[end]))
        end
        push!(
            records,
            Dict("source" => source, "n" => n, "twist" => twist, "entries_checked" => length(h)),
        )
    end
    return records
end

function triad_checks()
    weight(x) = x[1] == x[2] == x[3] ? 4 : x[1] == x[2] || x[2] == x[3] ? 3 : 2
    records = Any[]
    for xx in Iterators.product(0:2, 0:2, 0:2)
        x = Tuple(xx)
        row = Dict{NTuple{3,Int},Int}()
        for j in 1:2
            k = j + 1
            if x[j] == x[k]
                for s in 0:2
                    s == x[j] && continue
                    y = ntuple(l -> l in (j, k) ? s : x[l], 3)
                    row[y] = get(row, y, 0) + 1
                end
            else
                y = ntuple(l -> l == j ? x[k] : l == k ? x[j] : x[l], 3)
                row[y] = get(row, y, 0) - 1
            end
        end
        bound = sum(abs(v) * weight(y) for (y, v) in row)
        @assert bound <= 3weight(x)
        push!(records, Dict("word" => x, "weight" => weight(x), "weighted_row_sum" => bound))
    end
    return records
end

function error_bound_checks()
    @assert big(98)^281 < factorial(big(281))
    @assert frac(98, 282) < frac(1, 2) && frac(1, big(2)^96) < frac(1, big(2)^55)
    c1, q1 = poisson_coefficients(1)
    c2, q2 = poisson_coefficients(2)
    @assert all(>=(0), c1) && all(>=(0), c2)
    records = Any[]
    for n in 4:12
        d = isqrt(big(3)^n) + 1
        error = 281 * (d + 1)
        @assert d^2 > 3^n && q1 + error < 2q1
        row_bound = abs(32 - 3n) + 6n
        @assert row_bound <= 76 && row_bound * 2q1 < big(2)^63
        @assert cld(row_bound * 2q1, 32) + q1 < big(2)^63
        x = frac(49, 4) * n * (frac(3, 2) - frac(700741, 500000))
        @assert 0 <= x < 15 && x / (101 - x) < 1
        push!(
            records,
            Dict(
                "n" => n,
                "column_error_Q_units" => string(error + 1),
                "frobenius_error_Q_units" => string(d * (error + 1)),
                "int64_numerator_bound" => string(row_bound * 2q1),
            ),
        )
    end
    error2 = 257 * (2 * 243 + 1)
    @assert error2 == 125159
    @assert 32 * (q2 + error2) < big(2)^63 && 2 * (q2 + error2) < big(2)^63
    @assert big(27)^28 < big(32)^28
    for n in (6, 8, 10)
        x = frac(98518, 10^6) * frac(21n, 2)
        @assert 0 <= x < 12 && x / 102 < frac(1, 2)
    end
    # The largest possible term is <= 3*(2*(Q2+error2))^2. This also bounds
    # the checked 256-term Int128 chunks used by the unpacked implementation.
    @assert 256 * 3 * (2 * (q2 + error2))^2 < typemax(Int128)
    return records
end

function run_arithmetic(output)
    result = Dict(
        "status" => "reproduced",
        "source_integrity" => verify_source(),
        "matrices" => matrix_and_symmetry_checks(),
        "triad" => triad_checks(),
        "bounds" => error_bound_checks(),
        "completed_utc" => string(now(UTC)),
        "method" => "exact finite checks of physical operators, symmetries and error constants",
    )
    println("PASS physical operators, twisted symmetries, all 27 triad rows, and arithmetic bounds")
    return write_result(output, result)
end
