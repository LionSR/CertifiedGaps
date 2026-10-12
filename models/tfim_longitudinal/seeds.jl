# Write results/seeds.json: the exact seed data of all 30 certificates in one small file,
# so downstream users (e.g. a Lean statement of the criterion) need not decompress spectra.
include("certify.jl")
using .LongitudinalCertificate, CertifiedGaps.Results
const C=LongitudinalCertificate

"Exact seed rows (g, h, n₀, β₀, εₛ, εₜ, bound) read from the per-point certificates."
function seed_rows()
    rows=Dict{String,Any}[]
    for path in C.result_paths()
        d=C.read_point(path); c=d["certificate"]
        push!(rows,Dict{String,Any}("g"=>d["g"],"h"=>d["h"],"n0"=>c["n"],"beta0"=>c["beta"],
            "epsilon_s"=>c["epsilon_s"],"epsilon_t"=>c["epsilon_t"],
            "bound_exact_binary"=>c["bound_exact_binary"],"file"=>basename(path)))
    end
    sort!(rows;by=r->(C.rat(r["g"]),C.rat(r["h"])))
end

if abspath(PROGRAM_FILE)==@__FILE__
    write_result(joinpath(C.ROOT,"results","seeds.json"),Dict(
        "hamiltonian"=>"H_L = -sum_j Z_j Z_{j+1} - g sum_j X_j - h sum_j Z_j (Pauli, periodic)",
        "rp"=>"-h_bond = Z (x) Z with C = Z real; onsite o = -g X - h Z real symmetric",
        "epsilon_s"=>"upper bound on 1 - Z_{2n0}(beta0)/Z_{n0}(beta0)^2",
        "epsilon_t"=>"upper bound on 1 - Z_{2n0}(2 beta0)/Z_{2n0}(beta0)^2",
        "validity"=>"every even L >= n0 sites","points"=>seed_rows()))
end
