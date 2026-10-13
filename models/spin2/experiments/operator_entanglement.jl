# Numerical operator Schmidt spectra of the thermal amplitude A(β) = e^{−βH/2} on an OPEN spin-2
# chain, at the inverse temperatures of a squaring schedule β = 2^k. For each β the central-cut
# spectrum gives the bond dimension a compressed MPO needs for a given relative HS truncation
# error, which is the cost driver of the certified-squaring route (IDEAS.md, route H).
# Floating point, truncated at `maxdim`, open boundaries (edge spin-1 modes present at large β):
# an estimate, never a certificate.
#
#   julia --startup-file=no --project=models/spin2 models/spin2/experiments/operator_entanglement.jl
using ITensors, ITensorMPS, LinearAlgebra, JSON3, Dates
using Spin2Gap
using CertifiedGaps.Results: write_result

function operator_entanglement(;
    n = 32,
    maxdim = 200,
    cutoff = 1e-12,
    output = joinpath(@__DIR__, "operator-entanglement-n$(n)-d$(maxdim).json"),
)
    sites = siteinds("S=2", n; conserve_qns = true)
    pairs = bonds(n; periodic = false)
    A = MPO(sites, "Id")
    for j in 1:n
        A[j] /= sqrt(5)
    end
    # Schedule: step 1/32 to β = 2, then 1/8 to β = 8, then 1/2 to β = 64; targets are 2^k.
    schedule = [(1 / 32, 2.0), (1 / 8, 8.0), (1 / 2, 64.0)]
    targets = Set(2.0^k for k in -5:6)
    centre = n ÷ 2
    rows = Any[]
    β, lognorm, start = 0.0, 0.0, time()
    for (dβ, βend) in schedule
        halfgates = [exp(-(dβ / 4) * Spin2Gap.bond_operator(sites, i, j)) for (i, j) in pairs]
        gates = vcat(halfgates, reverse(halfgates))
        while β < βend - 1e-12
            A = apply(gates, A; cutoff, maxdim)
            g = norm(A)
            lognorm += log(g)
            A /= g
            β += dβ
            any(t -> isapprox(t, β; atol = 1e-9), targets) || continue
            A = orthogonalize(A, centre)
            left = uniqueinds(A[centre], A[centre+1])
            _, _, _, spec = svd(A[centre], left; cutoff = 0.0, maxdim = 10^6)
            p = sort(eigs(spec); rev = true)            # squared operator Schmidt values
            p ./= sum(p)
            tail = reverse(cumsum(reverse(p)))           # tail[D+1] = Σ_{i>D} p_i
            discarded = Dict(
                string(D) => (D < length(p) ? tail[D+1] : 0.0) for
                D in (8, 16, 32, 64, 128, 256, 512)
            )
            entropy = -sum(x -> x > 0 ? x * log(x) : 0.0, p)
            row = Dict(
                "beta" => β,
                "logZ_estimate" => n * log(5.0) + 2lognorm,
                "maxlinkdim" => maxlinkdim(A),
                "kept" => length(p),
                "operator_entropy" => entropy,
                "discarded_weight_at_D" => discarded,
                "seconds" => time() - start,
                "spectrum_top" => p[1:min(end, 40)],
            )
            push!(rows, row)
            println(
                "β=$β chi=$(maxlinkdim(A)) S_op=$(round(entropy; digits = 4)) ",
                "w(32)=$(discarded["32"]) w(64)=$(discarded["64"]) w(128)=$(discarded["128"]) ",
                "t=$(round(time() - start; digits = 1))s",
            )
            flush(stdout)
            write_result(
                output,
                Dict(
                    "status" => "numerical-operator-entanglement",
                    "certified" => false,
                    "n" => n,
                    "maxdim" => maxdim,
                    "cutoff" => cutoff,
                    "boundary" => "open",
                    "rows" => rows,
                    "created_utc" => string(now(UTC)),
                    "note" => "discarded weights at D >= maxlinkdim are understated",
                ),
            )
        end
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    BLAS.set_num_threads(2)
    operator_entanglement(;
        n = parse(Int, get(ENV, "OE_N", "32")),
        maxdim = parse(Int, get(ENV, "OE_MAXDIM", "200")),
    )
end
