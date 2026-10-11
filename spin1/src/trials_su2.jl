# Exact SU(2) multiplicity reduction. Every charge block is reconstructed by an
# invertible rational intertwiner, so correctness is a finite matrix identity.
import Nemo

function virtual_ladders(labels, A)
    plus, minus = zeros(BigInt, 26, 26), zeros(BigInt, 26, 26)
    lookup = Dict(label => i for (i, label) in enumerate(labels))
    for (j, (l, a, d)) in enumerate(labels)
        r = div(l + d, 2)
        d < l && (plus[lookup[(l, a, d + 2)], j] = l - r)
        d > -l && (minus[lookup[(l, a, d - 2)], j] = r)
    end
    @assert plus * minus - minus * plus == Matrix(Diagonal(big.([label[3] for label in labels])))
    for (generator, s, target, factor) in (
        (plus, -1, -1, 0),
        (plus, 0, -1, -1),
        (plus, 1, 0, -1),
        (minus, -1, 0, -2),
        (minus, 0, 1, -2),
        (minus, 1, 1, 0),
    )
        @assert generator * A[s] - A[s] * generator == factor * A[target]
    end
    return plus, minus
end

function su2_trial_blocks(labels, A, charge_blocks)
    plus, minus = virtual_ladders(labels, A)
    pairs = Dict(
        delta => [
            (i, j) for i in eachindex(labels) for
            j in eachindex(labels) if labels[i][3] - labels[j][3] == delta
        ] for delta in -10:2:10
    )
    lookup = Dict(
        delta => Dict(pair => i for (i, pair) in enumerate(pairs[delta])) for delta in keys(pairs)
    )
    raising, lowering = Dict{Int,Nemo.QQMatrix}(), Dict{Int,Nemo.QQMatrix}()
    for delta in -10:2:10
        if delta < 10
            R = zeros(BigInt, length(pairs[delta+2]), length(pairs[delta]))
            for (column, (a, b)) in enumerate(pairs[delta]), c in eachindex(labels)
                !iszero(plus[c, a]) && (R[lookup[delta+2][(c, b)], column] += 2plus[c, a])
                !iszero(minus[c, b]) && (R[lookup[delta+2][(a, c)], column] -= minus[c, b])
            end
            raising[delta] = Nemo.matrix(Nemo.QQ, R)
        end
        if delta > -10
            L = zeros(BigInt, length(pairs[delta-2]), length(pairs[delta]))
            for (column, (a, b)) in enumerate(pairs[delta]), c in eachindex(labels)
                !iszero(minus[c, a]) && (L[lookup[delta-2][(c, b)], column] += minus[c, a])
                !iszero(plus[c, b]) && (L[lookup[delta-2][(a, c)], column] -= 2plus[c, b])
            end
            lowering[delta] = Nemo.matrix(Nemo.QQ, L)
        end
    end
    T = Dict(b.delta => Nemo.matrix(Nemo.QQ, b.T) for b in charge_blocks)
    O = Dict(b.delta => Nemo.matrix(Nemo.QQ, b.O) for b in charge_blocks)
    highest = Dict{Int,Nemo.QQMatrix}()
    reduced = Any[]
    for spin in 0:5
        delta = 2spin
        multiplicity, N =
            spin == 5 ? (1, Nemo.identity_matrix(Nemo.QQ, 1)) : Nemo.nullspace(raising[delta])
        expected = [26, 56, 48, 24, 7, 1][spin+1]
        @assert multiplicity == expected
        @assert Nemo.rank(N) == multiplicity
        # RREF of N^t identifies independent rows; solve the square restriction.
        _, echelon = Nemo.rref(transpose(N))
        pivots =
            [findfirst(j -> !iszero(echelon[i, j]), 1:Nemo.ncols(echelon)) for i in 1:multiplicity]
        left_inverse = inv(N[pivots, 1:multiplicity])
        t = left_inverse * (T[delta]*N)[pivots, 1:multiplicity]
        o = left_inverse * (O[delta]*N)[pivots, 1:multiplicity]
        @assert T[delta] * N == N * t && O[delta] * N == N * o
        highest[spin] = N
        push!(reduced, (spin = spin, dimension = multiplicity, T = t, O = o))
    end
    # Explicitly reconstruct EVERY original block, not only the highest weights.
    # This verifies trace preservation for every power, without a sampled test.
    for delta in -10:2:10
        columns, timages, oimages = Nemo.QQMatrix[], Nemo.QQMatrix[], Nemo.QQMatrix[]
        for spin in div(abs(delta), 2):5
            N = highest[spin]
            for weight in (2spin):-2:(delta+2)
                N = lowering[weight] * N
            end
            push!(columns, N)
            push!(timages, N * reduced[spin+1].T)
            push!(oimages, N * reduced[spin+1].O)
        end
        W = reduce(hcat, columns)
        @assert size(W) == size(T[delta]) && Nemo.rank(W) == size(W, 1)
        @assert T[delta] * W == reduce(hcat, timages)
        @assert O[delta] * W == reduce(hcat, oimages)
    end
    @assert sum((2b.spin + 1) * b.dimension for b in reduced) == 676
    return reduced
end

function integral_scale(A::Nemo.QQMatrix)
    denominator = foldl(
        lcm,
        (BigInt(Nemo.denominator(A[i, j])) for i in 1:size(A, 1) for j in 1:size(A, 2));
        init = big(1),
    )
    return Nemo.matrix(Nemo.ZZ, A * denominator), denominator
end

function contract_trial_su2(reduced, n)
    irreps = Any[]
    for block in reduced
        t, dt = integral_scale(block.T)
        o, do_ = integral_scale(block.O)
        X = t^(n - 2)
        u_numerator = BigInt(Nemo.tr(X * o))
        v_numerator = BigInt(Nemo.tr(X * t^2))
        du, dv = dt^(n - 2) * do_, dt^n
        @assert iszero(rem(u_numerator, du)) && iszero(rem(v_numerator, dv))
        push!(
            irreps,
            (
                spin = block.spin,
                U = div(u_numerator, du),
                V = div(v_numerator, dv),
                dimension = block.dimension,
            ),
        )
    end
    # A spin-J copy contributes once to each weight -J,...,J.
    records = Any[]
    for delta in -10:2:10
        selected = filter(b -> 2b.spin >= abs(delta), irreps)
        push!(
            records,
            Dict(
                "delta" => delta,
                "dimension" => sum(b.dimension for b in selected),
                "U" => string(sum(b.U for b in selected)),
                "V" => string(sum(b.V for b in selected)),
            ),
        )
    end
    u = sum((2b.spin + 1) * b.U for b in irreps)
    v = sum((2b.spin + 1) * b.V for b in irreps)
    return u, v, records
end

function run_trials(output; backend = :su2)
    backend == :reference && return run_trials_reference(output)
    backend == :su2 || throw(ArgumentError("trial backend must be :su2 or :reference"))
    provenance = verify_source()
    started = time()
    labels, A = trial_matrices()
    blocks = trial_blocks(labels, A; native = true)
    physical = physical_trial_check(A, blocks)
    reduced = su2_trial_blocks(labels, A, blocks)
    reduction_seconds = time() - started
    records = Any[]
    for n in (72, 120)
        u, v, block_records = contract_trial_su2(reduced, n)
        expected = read_exact_json(joinpath(SOURCE, "evidence/trials/certificate$n.json"))
        @assert u == integer(expected[n == 72 ? "one_bond_energy_numerator" : "U"])
        @assert v == integer(expected[n == 72 ? "norm" : "V"])
        margin = -700741v - 500000u
        @assert v > 0 && margin > 0
        @assert margin == integer(
            expected[n == 72 ? "positive_integer_margin" :
                     "positive_margin_minus700741V_minus500000U"],
        )
        if n == 120
            for (actual, reference) in zip(block_records, expected["blocks"]),
                key in ("delta", "dimension", "U", "V")

                @assert integer(actual[key]) == integer(reference[key])
            end
        end
        push!(
            records,
            Dict(
                "length" => n,
                "U" => string(u),
                "V" => string(v),
                "positive_margin" => string(margin),
                "floor_1e12_U_over_V" => string(fld(big(10)^12 * u, v)),
                "matches_upstream_exactly" => true,
                "blocks" => block_records,
            ),
        )
        println(
            "PASS SU(2) exact trial N=$n; full U,V and charge blocks agree; elapsed=$(round(time()-started; digits=3))s",
        )
        flush(stdout)
    end
    return write_result(
        output,
        Dict(
            "status" => "reproduced",
            "source_integrity" => provenance,
            "method" => "exact SU(2) multiplicity blocks with verified rational similarity and FLINT integer powers",
            "multiplicity_dimensions" => [b.dimension for b in reduced],
            "physical_four_site_check" => physical,
            "reduction_seconds" => reduction_seconds,
            "wall_seconds" => time() - started,
            "julia" => string(VERSION),
            "Nemo" => string(pkgversion(Nemo)),
            "completed_utc" => string(now(UTC)),
            "certificates" => records,
        ),
    )
end
