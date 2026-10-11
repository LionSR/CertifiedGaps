import Spin1Certificates as SC
basis = SC.magnetization_basis(Val(12), 0)
reps, mult = SC.representatives(basis, 1, 0)
rows = SC.thermal_rows(basis, 1, 0)
coefficients, _ = SC.poisson_coefficients(1)

function partial_contraction(batch, rows, reps, mult, coefficients; columns = 64)
    started = time()
    total = big(0)
    for first in 1:batch:columns
        last = min(first + batch - 1, columns)
        rr = reps[first:last]
        u = zeros(Int64, length(rr), length(rows))
        y = similar(u)
        for q in Iterators.reverse(coefficients)
            SC.multiply_real!(y, rows, u, 5)
            for (lane, row) in enumerate(rr)
                y[lane, row] += q
            end
            u, y = y, u
        end
        for (lane, multiplicity) in enumerate(mult[first:last])
            total += multiplicity * SC.squared_norm(u, u, lane, false)
        end
    end
    return total, time() - started
end

function main()
    partial_contraction(4, rows, reps, mult, coefficients; columns = 4)
    records = Any[]
    reference = nothing
    for batch in (4, 8, 16, 32)
        GC.gc()
        total, seconds = partial_contraction(batch, rows, reps, mult, coefficients)
        isnothing(reference) ? (reference = total) : (@assert total == reference)
        push!(
            records,
            Dict(
                "batch" => batch,
                "seconds" => seconds,
                "two_buffer_bytes" => 2 * 8 * batch * length(rows),
                "partial_total" => string(total),
            ),
        )
        println("n=12 Sz=0 first 64 columns: batch=$batch seconds=$seconds")
        flush(stdout)
    end
    SC.write_result(
        joinpath(@__DIR__, "../results/batch-benchmark.json"),
        Dict(
            "status" => "passed",
            "scope" => "Same 64 columns and all 192 Horner steps of the n=12 Sz=0 sector; microbenchmark, not a complete certificate",
            "comparisons" => records,
        ),
    )
end
main()
