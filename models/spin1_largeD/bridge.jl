isdefined(@__MODULE__, :TransferCertificate) || include("src/TransferCertificate.jl")
using CertifiedGaps.Results

"Fresh full-spectrum enclosures for the two small rings omitted by the n0=12 seed."
function bridge_certificate()
    tc = TransferCertificate
    started = time()
    rows = Any[]
    for L in (8, 10)
        s = tc.compute_spectrum(L, 4//1; rigorous = true)
        gap = tc.finite_gap((; E = s["energies"], R = s["radii"], L, D = 4//1), 4//1, 4//1)
        tc.inf(gap) > 0 || error("finite gap not separated")
        # A coarse exact enclosure tolerates harmless LAPACK/BLAS replay differences.
        grid = big(10)^6
        lower = floor(BigInt, tc.Q(tc.inf(gap)) * grid) // grid
        higher = ceil(BigInt, tc.Q(tc.sup(gap)) * grid) // grid
        push!(rows, Dict("L" => string(L), "D" => "4//1",
            "gap_interval" => string.([lower, higher]), "dimension" => string(3^L),
            "blocks" => string(length(s["sectors"])), "seconds" => s["seconds"]))
    end
    hashes = tc.proof_hashes()
    return Dict("rings" => rows, "wall_seconds" => time()-started,
        "source_hashes" => hashes,
        "run_metadata" => run_metadata(Dict("D" => "4//1");
            source_sha256 = hashes, packages = (tc.CertifiedGaps, tc.IntervalArithmetic),
            status = "fresh complete spectra; checker repeats these small eigensolves",
            certified = true))
end

"Recompute the small-ring bridge and require it to exceed the larger-seed bound."
function check_bridge(data, threshold)
    tc = TransferCertificate
    Dict(data["source_hashes"]) == tc.proof_hashes() || error("bridge source mismatch")
    fresh = bridge_certificate()
    length(data["rings"]) == 2 || error("incomplete finite bridge")
    for (stored, recomputed) in zip(data["rings"], fresh["rings"])
        stored["L"] == recomputed["L"] || error("wrong bridge length")
        stored["D"] == recomputed["D"] || error("wrong bridge parameter")
        stored["dimension"] == recomputed["dimension"] || error("wrong bridge dimension")
        stored["blocks"] == recomputed["blocks"] || error("wrong bridge block count")
        sl, sh = parse.(tc.Q, stored["gap_interval"])
        rl, rh = parse.(tc.Q, recomputed["gap_interval"])
        sl <= rl <= rh <= sh || error("bridge did not replay")
        parse(tc.Q, stored["gap_interval"][1]) >= threshold || error("bridge gap too small")
    end
    a = tc.reference_spectrum(4//1, 6)
    tc.Q(tc.inf(tc.finite_gap(a, 4//1, 4//1))) >= threshold || error("L=6 gap too small")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    write_result(joinpath(@__DIR__, "results", "finite_bridge.json"), bridge_certificate())
end
