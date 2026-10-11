# Compare candidate generators on exactly the same saved spin-2 tensor input.
# These timings do not constitute a thermal trace or uniform-gap certificate.
import CertifiedGaps
using CertifiedGaps.Results: write_result
using LinearAlgebra, JSON3, SHA
import Nemo
const MP = CertifiedGaps.CertifiedMPO.Products
const MS = MP.MS
const MN = MP.MN
const CB = MS.CB
const F = MP.F
BLAS.set_num_threads(1)

function benchmark_charge_blocks(source, output; repetitions = 3, maxdim = 34)
    @assert abspath(source) != abspath(output)
    source_text = read(source, String)
    data = JSON3.read(source_text)
    A = MN.decode_cores(data["candidate"])
    physical, weights = collect(-2:2), F[1, 4, 6, 4, 1]
    charges = CB.infer_charges(A, physical)
    identity = [reshape(Matrix{F}(I, 5, 5), 1, 5, 5, 1) for _ in A]
    RR = Nemo.ArbField(192)
    # Warm both paths on a small prefix-shaped example before timing.
    tiny = identity[1:2]
    for symmetry in (nothing, physical)
        MP.product_candidate(
            tiny,
            tiny;
            physical_charges = symmetry,
            candidate_arithmetic = :interval,
            maxdim,
        )
    end
    rows = Any[]
    for (name, symmetry) in (("dense", nothing), ("charge_blocks", physical))
        GC.gc()
        println("Generating $name candidate, N=$(length(A)), maxdim=$maxdim")
        flush(stdout)
        elapsed = @timed MP.product_candidate(
            A,
            identity;
            physical_charges = symmetry,
            weights,
            candidate_arithmetic = :interval,
            maxdim,
            bits = 64,
        )
        proof = elapsed.value
        println("$name generation=$(elapsed.time)s; replaying with the general checker")
        flush(stdout)
        replay_seconds =
            @elapsed @assert MP.verify_product(A, identity, proof; residual_arithmetic = :interval)
        witness = splitext(output)[1]*"-$name-witness.json"
        open(witness, "w") do io
            JSON3.write(io, MP.encode(A, identity, proof))
        end
        nonzero = sum(count(!iszero, c) for c in proof.cores)
        entries = sum(length, proof.cores)
        push!(
            rows,
            Dict(
                "algorithm" => name,
                "generation_seconds" => elapsed.time,
                "julia_allocated_bytes" => elapsed.bytes,
                "replay_seconds" => replay_seconds,
                "hs_error_upper" => string(proof.error),
                "hs_error_display" => Float64(proof.error),
                "output_nonzeros" => nonzero,
                "output_entries" => entries,
                "witness_file" => basename(witness),
                "witness_bytes" => filesize(witness),
                "witness_sha256" => bytes2hex(sha256(read(witness))),
                "general_interval_replay_passed" => true,
            ),
        )
        symmetry === nothing || CB.infer_charges(proof.cores, physical)
        println(
            "$name: error=$(Float64(proof.error)), density=$(nonzero/entries), replay=$(replay_seconds)s",
        )
        flush(stdout)
    end
    # Representative Gram multiplication: a saved physical core with a
    # synthetic charge-compatible bridge. Explicitly separated from full-run timings.
    site = min(4, length(A)-1)
    first = MN.adjoint_cores(A; weights)
    qa = [-q for q in charges]
    prefix = charges[site]
    bridge = [
        p == a+b ? F(mod(i+3j+7k, 11)-5, 16) : F(0) for (i, p) in enumerate(prefix),
        (j, a) in enumerate(qa[site]),
        (k, b) in enumerate(charges[site])
    ]
    dense() = MS.product_bridge(bridge, first[site], A[site]; field = RR)
    blocked() = CB.product_bridge(
        bridge,
        first[site],
        A[site],
        prefix,
        qa[site],
        qa[site+1],
        charges[site],
        charges[site+1],
        physical;
        field = RR,
    )
    reference = dense()
    trial = blocked()
    @assert all(Nemo.overlaps(a, b) for (a, b) in zip(reference, trial))
    kernel_rows = Any[]
    for (name, run) in (("dense", dense), ("charge_blocks", blocked))
        timings, allocations = Float64[], Int[]
        for _ in 1:repetitions
            GC.gc()
            elapsed = @timed run()
            push!(timings, elapsed.time)
            push!(allocations, elapsed.bytes)
        end
        push!(
            kernel_rows,
            Dict("algorithm" => name, "seconds" => timings, "julia_allocated_bytes" => allocations),
        )
    end
    result = Dict(
        "status" => "charge-block-optimization-benchmark",
        "source_file" => abspath(source),
        "source_sha256" => bytes2hex(sha256(source_text)),
        "input_nonzero_density" => sum(count(!iszero, c) for c in A)/sum(length, A),
        "length" => length(A),
        "maxdim" => maxdim,
        "bits" => 64,
        "precision_bits" => 192,
        "blas_threads" => BLAS.get_num_threads(),
        "julia_version" => string(VERSION),
        "nemo_version" => string(Base.pkgversion(Nemo)),
        "generator_runs" => rows,
        "gram_kernel_runs" => kernel_rows,
        "gram_kernel_bridge" => "synthetic charge-compatible bridge; saved physical cores",
        "kernel_interval_overlap_passed" => true,
        "uniform_gap_certified" => false,
        "thermal_trace_certified" => false,
    )
    write_result(output, result)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    @assert length(ARGS) == 2
    benchmark_charge_blocks(ARGS...)
end
