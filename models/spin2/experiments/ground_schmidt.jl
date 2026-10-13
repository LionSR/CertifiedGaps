# Numerical ground-state Schmidt spectrum of an OPEN spin-2 chain (DMRG). At β ≫ 1/Δ the
# operator Schmidt spectrum of e^{−βH/2} approaches the products λ_a λ_b of these values, which
# sets the large-β bond dimension of the certified-squaring route. Not a certificate.
#
#   julia --startup-file=no --project=models/spin2 models/spin2/experiments/ground_schmidt.jl
using ITensors, ITensorMPS, LinearAlgebra, JSON3, Dates
using Spin2Gap
using CertifiedGaps.Results: write_result

function ground_schmidt(;
    n = 64,
    maxdims = [50, 100, 200, 300, 400, 400, 400, 400],
    output = joinpath(@__DIR__, "ground-schmidt-n$(n).json"),
)
    sites = siteinds("S=2", n; conserve_qns = true)
    H = hamiltonian(sites; periodic = false)
    ψ0 = Spin2Gap.initial_state(sites, 0)                   # Sz = 0 Néel start, as in dmrg.jl
    nsweeps = length(maxdims)
    cutoff = 1e-12
    start = time()
    E, ψ = dmrg(H, ψ0; nsweeps, maxdim = maxdims, cutoff, outputlevel = 1)
    centre = n ÷ 2
    ψ = orthogonalize(ψ, centre)
    left = uniqueinds(ψ[centre], ψ[centre+1])
    _, _, _, spec = svd(ψ[centre], left; cutoff = 0.0, maxdim = 10^6)
    p = sort(eigs(spec); rev = true)
    p ./= sum(p)
    tail = reverse(cumsum(reverse(p)))
    discarded = Dict(
        string(D) => (D < length(p) ? tail[D+1] : 0.0) for D in (16, 32, 50, 100, 200, 300, 400)
    )
    result = Dict(
        "status" => "numerical-ground-schmidt",
        "certified" => false,
        "n" => n,
        "boundary" => "open",
        "energy" => E,
        "energy_per_bond" => E / (n - 1),
        "kept" => length(p),
        "entropy" => -sum(x -> x > 0 ? x * log(x) : 0.0, p),
        "discarded_weight_at_D" => discarded,
        "spectrum" => p,
        "maxdims" => maxdims,
        "seconds" => time() - start,
        "created_utc" => string(now(UTC)),
    )
    write_result(output, result)
    println("E=$E per bond $(E/(n-1)) entropy=$(result["entropy"]) discarded=$discarded")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    BLAS.set_num_threads(2)
    ground_schmidt(; n = parse(Int, get(ENV, "GS_N", "64")))
end
