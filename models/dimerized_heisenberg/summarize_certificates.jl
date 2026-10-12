# Writes results/certificate.json: the seed data of the three verified certificates, as plain
# exact rationals, without the full spectra. It only copies fields from the compressed
# verified_delta_*.json.gz files; check_verified.jl confirms the copy agrees with them.
using CertifiedGaps, CertifiedGaps.Results, JSON3
include("provenance.jl")

"Collect the seed data and RP description of every stored verified certificate."
function certificate_summary()
    dir = joinpath(@__DIR__, "results")
    paths = sort(filter(p -> startswith(basename(p), "verified_delta_") &&
        endswith(p, ".json.gz"), readdir(dir; join=true)))
    rows = Any[]
    for path in paths
        data = JSON3.read(read_json_text(path))
        push!(rows, Dict(
            "delta" => String(data.delta), "n0_cells" => Int(data.n0_cells),
            "n0_spins" => 2 * Int(data.n0_cells), "beta0" => String(data.beta0),
            "epsilon_s_upper" => String(data.epsilon_s_upper),
            "epsilon_t_upper" => String(data.epsilon_t_upper),
            "log_partition_enclosures" => Dict(String(k) => String.(v)
                for (k, v) in pairs(data.log_partition_enclosures)),
            "recurrence_bound_exact_binary" => String(data.recurrence_bound_exact_binary),
            "validity_range" => String(data.validity_range),
            "source_file" => "results/" * basename(path),
            "source_file_sha256" => sha256_file(path)))
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    root = repository_root()
    hashes = source_hashes(root; dirs=(), files=(
        "models/dimerized_heisenberg/summarize_certificates.jl",
        "models/dimerized_heisenberg/provenance.jl"))
    write_result(joinpath(@__DIR__, "results", "certificate.json"), Dict(
        "metadata" => run_metadata(Dict("model" => "dimerized_heisenberg",
                "sites_per_cell" => 2); source_sha256=hashes, packages=(CertifiedGaps,),
            status="summary_of_verified_certificates", certified=true),
        "model" => "H = sum over a ring of 2n spins of (1+delta) S.S inside each two-spin cell " *
            "plus (1-delta) S.S between neighbouring cells",
        "criterion_lengths" => "n counts two-spin cells; P_s(n0,beta0) = Z_{2n0}/Z_{n0}^2, " *
            "P_tau(2n0,beta0) = Z_{2n0}(2beta0)/Z_{2n0}(beta0)^2",
        "rp_form" => "Scope-note bond form with an internal reflection: rotate the second spin of " *
            "every cell by J = [0 1; -1 0]; then -h = (1-delta) sum_{C in (Sx, iSy, Sz) on the " *
            "second spin of the left cell} C (x) Theta(C), Theta(X) = F conj(X) F with F the " *
            "swap of the two spins of a cell, and o = (1+delta) J-rotated S1.S2 is real with " *
            "Theta(o) = o. In the cell basis M/sqrt(2) (F = M M^T/2, columns |00>+|11>, " *
            "|00>-|11>, |01>+|10>, i(|01>-|10>)) this is -h = sum C' (x) conj(C') with complex " *
            "C' and real diagonal o. No common cell basis gives a swap-symmetric sum C (x) C. " *
            "Exact checks: rp.jl.",
        "certificates" => certificate_summary()))
end
