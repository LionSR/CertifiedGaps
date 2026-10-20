# Benchmark certified Arb multiplication on shapes used by `product_bridge`.
# The matrices are synthetic; these timings are not thermal-certificate inputs.
# Both routines retain rigorous error bounds, as documented in
# https://flintlib.org/doc/arb_mat.html. The block routine uses integer products.
import Nemo
using JSON3
using SHA
using CertifiedGaps.Results: write_result

function block_product(A::Nemo.ArbMatrix, B::Nemo.ArbMatrix)
    @assert Nemo.ncols(A) == Nemo.nrows(B)
    @assert Nemo.base_ring(A) === Nemo.base_ring(B)
    C = similar(A, Nemo.nrows(A), Nemo.ncols(B))
    ccall(
        (:arb_mat_mul_block, Nemo.libflint),
        Cvoid,
        (Ref{Nemo.ArbMatrix}, Ref{Nemo.ArbMatrix}, Ref{Nemo.ArbMatrix}, Int),
        C,
        A,
        B,
        precision(Nemo.base_ring(A)),
    )
    return C
end

function synthetic_product_inputs(field, m, k, n)
    # Include signed entries, variation in scale, and non-dyadic inputs whose
    # Arb conversion has a nonzero radius. All source entries are exact QQ.
    entry(i, j, salt) = Nemo.QQ(mod(17i+31j+salt, 101)-50, 3*big(2)^mod(7i+11j+salt, 24))
    A = Nemo.matrix(field, [entry(i, j, 0) for i in 1:m, j in 1:k])
    B = Nemo.matrix(field, [entry(i, j, 5) for i in 1:k, j in 1:n])
    return A, B
end

function benchmark_products(output; repetitions = 2, precision_bits = 192)
    @assert repetitions >= 1
    RR = Nemo.ArbField(precision_bits)
    exact_A, exact_B = synthetic_product_inputs(Nemo.QQ, 12, 7, 9)
    exact_C = exact_A*exact_B
    small_A, small_B = Nemo.matrix(RR, exact_A), Nemo.matrix(RR, exact_B)
    @assert Nemo.contains(small_A*small_B, exact_C)
    @assert Nemo.contains(block_product(small_A, small_B), exact_C)
    rows = Any[]
    for shape in ((512, 16, 400), (2560, 80, 80), (1024, 32, 800), (5120, 160, 160))
        m, k, n = shape
        A, B = synthetic_product_inputs(RR, m, k, n)
        baseline = A*B
        @assert Nemo.overlaps(block_product(A, B), baseline)
        for (name, multiply) in (("default", *), ("block", block_product))
            times = Float64[]
            for _ in 1:repetitions
                GC.gc()
                elapsed = @elapsed C = multiply(A, B)
                @assert Nemo.overlaps(C, baseline)
                push!(times, elapsed)
            end
            row = Dict("shape"=>[m, k, n], "algorithm"=>name, "seconds"=>times)
            push!(rows, row)
            println("shape=$shape algorithm=$name seconds=$times")
            flush(stdout)
        end
    end
    result = Dict(
        "status"=>"synthetic-certified-kernel-benchmark",
        "uniform_gap_certified"=>false,
        "precision_bits"=>precision_bits,
        "repetitions"=>repetitions,
        "small_exact_containment_passed"=>true,
        "large_overlap_checks_passed"=>true,
        "source_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "julia_version"=>string(VERSION),
        "nemo_version"=>string(Base.pkgversion(Nemo)),
        "rows"=>rows,
    )
    write_result(output, result)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 1
    benchmark_products(only(ARGS))
end
