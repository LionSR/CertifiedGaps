# Integer MPS recipe and transfer contractions from trials/original_trial.py.
function trial_matrices()
    multiplets = [(l, a) for (l, k) in ((1, 4), (3, 3), (5, 1)) for a in 0:(k-1)]
    labels = [(l, a, d) for (l, a) in multiplets for d in (-l):2:l]
    diagonal = Dict(zip(multiplets, [-289, -84, -9, 33, -22, 27, 71, -11]))
    pairs = [(x, y) for x in multiplets for y in multiplets if y[1] - x[1] == 2]
    off = Dict(zip(pairs, [-11, -35, -28, -23, 16, -49, -6, -25, -20, -34, 9, -32, 1, 2, -16]))
    lookup = Dict(label => i for (i, label) in enumerate(labels))
    A = Dict(s => zeros(BigInt, length(labels), length(labels)) for s in -1:1)
    for (i, (l, a, d)) in enumerate(labels), (k, b) in multiplets
        abs(k - l) <= 2 || continue
        p =
            k == l ? (a == b ? diagonal[(l, a)] : 0) :
            k > l ? off[((l, a), (k, b))] : -off[((k, b), (l, a))]
        r = div(d + l, 2)
        for s in -1:1
            abs(d + 2s) <= k || continue
            factors =
                k > l ?
                (2 * (l + 1 - r) * (l + 2 - r), 2 * (r + 1) * (l - r + 1), (r + 1) * (r + 2)) :
                k == l ? (2 * (l - r + 1), d, -(r + 1)) : (2, -2, 1)
            A[s][i, lookup[(k, b, d + 2s)]] = p * factors[s+2]
        end
    end
    expected = read_exact_json(joinpath(SOURCE, "evidence/trials/matrix-data.json"))
    @assert labels == [Tuple(Int.(integer.(label))) for label in expected["labels"]]
    for s in -1:1, i in eachindex(labels), j in eachindex(labels)
        @assert A[s][i, j] == integer(expected["matrices"][string(s)][i][j])
    end
    return labels, A
end

function trial_blocks(labels, A; native = false)
    if native
        amax = maximum(maximum(abs, a) for a in values(A))
        bmax = length(labels) * amax^2
        @assert 12bmax^2 < typemax(Int64)
        A = Dict(s => Int64.(a) for (s, a) in A)
    end
    element_type = eltype(first(values(A)))
    B = Dict((s, t) => A[s] * A[t] for s in -1:1 for t in -1:1)
    blocks = Any[]
    for delta in -10:2:10
        pairs = [
            (i, j) for i in eachindex(labels) for
            j in eachindex(labels) if labels[i][3] - labels[j][3] == delta
        ]
        T, O = zeros(element_type, length(pairs), length(pairs)),
        zeros(element_type, length(pairs), length(pairs))
        for (row, (a, b)) in enumerate(pairs), (column, (c, d)) in enumerate(pairs)
            T[row, column] = sum(A[s][a, c] * A[s][b, d] for s in -1:1)
            O[row, column] =
                sum(s * t * B[(s, t)][a, c] * B[(s, t)][b, d] for s in -1:1 for t in -1:1)
            for s in -1:0, t in 0:1
                O[row, column] +=
                    B[(s, t)][a, c] * B[(s + 1, t - 1)][b, d] +
                    B[(s + 1, t - 1)][a, c] * B[(s, t)][b, d]
            end
        end
        push!(blocks, (delta = delta, T = BigInt.(T), O = BigInt.(O)))
    end
    @assert [size(b.T, 1) for b in blocks] == [1, 8, 32, 80, 136, 162, 136, 80, 32, 8, 1]
    return blocks
end

function contract_trial(blocks, n)
    numerator, norm = big(0), big(0)
    records = Any[]
    for block in blocks
        X = block.T^(n - 2)
        u = sum(X .* transpose(block.O))
        v = sum(X .* transpose(block.T * block.T))
        numerator += u
        norm += v
        push!(
            records,
            Dict(
                "delta" => block.delta,
                "dimension" => size(block.T, 1),
                "U" => string(u),
                "V" => string(v),
            ),
        )
    end
    return numerator, norm, records
end

function physical_trial_check(A, blocks)
    words = [Tuple(w) for w in Iterators.product(ntuple(_ -> -1:1, 4)...)]
    amplitudes = Dict(word => tr(prod(A[s] for s in word)) for word in words)
    norm = sum(a^2 for a in values(amplitudes))
    energy_numerator = big(0)
    for (word, amp) in amplitudes
        sum(word) != 0 && @assert iszero(amp)
        @assert amp == amplitudes[(word[2], word[3], word[4], word[1])]
        for j in 1:4
            k = mod1(j + 1, 4)
            energy_numerator += word[j] * word[k] * amp^2
            for delta in (-1, 1)
                abs(word[j] - delta) <= 1 && abs(word[k] + delta) <= 1 || continue
                target = ntuple(l -> word[l] + (l == j ? -delta : l == k ? delta : 0), 4)
                energy_numerator += amp * amplitudes[target]
            end
        end
    end
    u, v, _ = contract_trial(blocks, 4)
    @assert norm == v == parse(BigInt, "22975205960525898253776")
    @assert energy_numerator == 4u == parse(BigInt, "-137851234864097615402880")
    return Dict(
        "length" => 4,
        "norm" => string(norm),
        "energy_numerator" => string(energy_numerator),
        "physical_and_transfer_agree" => true,
    )
end

function run_trials_reference(output)
    provenance = verify_source()
    start = time()
    labels, A = trial_matrices()
    blocks = trial_blocks(labels, A)
    physical = physical_trial_check(A, blocks)
    records = Any[]
    for n in (72, 120)
        u, v, block_records = contract_trial(blocks, n)
        margin = -700741v - 500000u
        @assert v > 0 && margin > 0
        expected = read_exact_json(joinpath(SOURCE, "evidence/trials/certificate$n.json"))
        @assert u == integer(expected[n == 72 ? "one_bond_energy_numerator" : "U"])
        @assert v == integer(expected[n == 72 ? "norm" : "V"])
        @assert margin == integer(
            expected[n == 72 ? "positive_integer_margin" :
                     "positive_margin_minus700741V_minus500000U"],
        )
        floor_energy = fld(big(10)^12 * u, v)
        @assert floor_energy == (n == 72 ? -1401482112601 : -1401482105174)
        if n == 120
            for (actual, reference) in zip(block_records, expected["blocks"])
                for key in ("delta", "dimension", "U", "V")
                    @assert integer(actual[key]) == integer(reference[key])
                end
            end
        end
        push!(
            records,
            Dict(
                "length" => n,
                "U" => string(u),
                "V" => string(v),
                "positive_margin" => string(margin),
                "floor_1e12_U_over_V" => string(floor_energy),
                "matches_upstream_exactly" => true,
                "blocks" => block_records,
            ),
        )
        println(
            "PASS integer MPS trial N=$n floor(10^12 U/V)=$floor_energy elapsed=$(round(time()-start; digits=1))s",
        )
        flush(stdout)
    end
    return write_result(
        output,
        Dict(
            "status" => "reproduced",
            "source_integrity" => provenance,
            "method" => "fresh BigInt MPS contractions from the polynomial entry recipe",
            "julia" => string(VERSION),
            "wall_seconds" => time() - start,
            "completed_utc" => string(now(UTC)),
            "physical_four_site_check" => physical,
            "certificates" => records,
        ),
    )
end
