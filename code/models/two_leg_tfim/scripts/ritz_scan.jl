# Floating-point diagnostic (not certified) of the coupled half-length seed at n₀ = 6 rungs:
#
#     εₛ = 1 − Z₁₂(β)/Z₆(β)²,
#     εₜ = 1 − [Z₁₂(2β)/Z₆(2β)²] · [Z₆(2β)/Z₆(β)²]²  ≥ 1 − P_τ(12, β)
#
# (coupling identity P_τ(2n, β) ≥ Pₛ(n, 2β) P_τ(n, β)², AGENTS.md). Z₆ is summed over the
# complete spectrum of the 6-rung ring (dense sectors); Z₁₂ is only bounded below, by the
# lowest Lanczos values of every sector of the 12-rung ring (24 spins, 28 sectors of about
# 350 000 states, scripts/sectors.jl): K₊ values in the flip-even sectors (ground state and
# two-quasiparticle states), K₋ in the flip-odd ones (one- and three-quasiparticle states).
# Omitting states only lowers Z₁₂, so truncating to the
# k lowest values per sector shows how many states a Ritz certificate needs.
#
#     julia --project=code/julia code/models/two_leg_tfim/scripts/ritz_scan.jl 5//2 [K₊ K₋]
#
# Writes results/ritz_spectra_g<g>.json (per-sector Lanczos values) and results/ritz_scan_g<g>.json.
const ROOT = dirname(@__DIR__)
include(joinpath(@__DIR__, "sectors.jl"))
using CertifiedGaps, CertifiedGaps.Results
using KrylovKit
BLAS.set_num_threads(1)

g = parse(Rational{Int}, get(ARGS, 1, "5//2"))
K = parse(Int, get(ARGS, 2, "14"))
Kodd = parse(Int, get(ARGS, 3, "6"))
tag = "$(numerator(g))_$(denominator(g))"
spectra_path = joinpath(ROOT, "results", "ritz_spectra_g$tag.json")

"Complete spectrum of the ring of L rungs as (values, multiplicities), dense per sector."
function full_spectrum(L, g)
    G = LadderGroup(L); tab = orbit_table(G); hops = hopping_table(G, tab)
    E, w = Float64[], Int[]
    for (m, s, f, mult) in sectors(L)
        A = sector_matrix(G, tab, hops, Float64(g), m, s, f)
        size(A, 1) == 0 && continue
        ev = eigvals(Hermitian(Matrix(A)))
        append!(E, ev); append!(w, fill(mult, length(ev)))
    end
    sum(w) == 4^L || error("sectors cover $(sum(w)) of $(4^L) states")
    return E, w
end

"The K (flip-even) or Kodd (flip-odd) lowest Lanczos values of every sector of the L-rung ring."
function low_spectrum(L, g, K, Kodd)
    G = LadderGroup(L); tab = orbit_table(G); hops = hopping_table(G, tab)
    out = []
    for (m, s, f, mult) in sectors(L)
        A = sector_matrix(G, tab, hops, Float64(g), m, s, f)
        x0 = ones(eltype(A), size(A, 1)) .+ 0.1 .* rand(eltype(A), size(A, 1))
        k = f == 1 ? K : Kodd
        seconds = @elapsed vals, _, info = eigsolve(A, x0, k, :SR; ishermitian = true,
            krylovdim = 2k + 10, tol = 1e-7, maxiter = 200)
        vals = sort(real.(vals))[1:min(k, end)]
        push!(out, (; m, s, f, multiplicity = mult, dimension = size(A, 1), values = vals,
            converged = info.converged, seconds))
        println("L=$L g=$g sector ($m, $s, $f) dim $(size(A, 1)) E₁ = $(vals[1]) ",
            "converged $(info.converged) $(round(seconds; digits = 1)) s"); flush(stdout)
    end
    return out
end

E6, w6 = full_spectrum(6, g)
low = if isfile(spectra_path) && !("--fresh" in ARGS)
    [(; multiplicity = Int(r.multiplicity), values = Vector{Float64}(r.values))
     for r in read_result(spectra_path; exact = false).sectors]
else
    s = low_spectrum(12, g, K, Kodd)
    write_result(spectra_path, (; status = "numerical, not certified (Lanczos)", g = string(g),
        K, Kodd, sectors = s); compact = true)
    s
end

# Shifted log-sum-exp so that large β does not underflow; e₀ is a common reference energy.
logZ(E, w, β, e0) = log(sum(w[i] * exp(-β * (E[i] - e0)) for i in eachindex(E)))
function logZ12(β, k, e0)
    E = reduce(vcat, [r.values[1:min(k, end)] for r in low])
    w = reduce(vcat, [fill(r.multiplicity, min(k, length(r.values))) for r in low])
    return logZ(E, w, β, e0)
end

function defects(β, k)
    e0 = minimum(E6) * 2                         # E₀(12) ≈ 2 E₀(6) keeps exponents moderate
    lz6, lz6_2 = logZ(E6, w6, β, e0 / 2), logZ(E6, w6, 2β, e0 / 2)
    lz12, lz12_2 = logZ12(β, k, e0), logZ12(2β, k, e0)
    εs = -expm1(lz12 - 2lz6)
    εt = -expm1((lz12_2 - 2lz6_2) + 2(lz6_2 - 2lz6))
    return εs, εt
end
bound(εs, εt, β) = (0 <= εs < 1 && 0 <= εt < 1) ? certified_gap_bound(εs, εt, β).bound :
                   nothing

βs = 1/8:1/8:16
ks = unique([1, 2, 3, 4, 6, 8, 12, K])
rows = []
for k in ks, β in βs
    εs, εt = defects(β, k)
    push!(rows, (; k, beta = β, eps_s = εs, eps_t = εt, bound = bound(εs, εt, β)))
end
best = Dict(k => let ok = filter(r -> r.k == k && r.bound !== nothing, rows)
    isempty(ok) ? nothing : argmax(r -> r.bound, ok)
end for k in ks)
gap12 = let v = sort(reduce(vcat, [r.values for r in low])); v[2] - v[1] end
for k in ks
    println("g = $g, k = $k states/sector: best ", best[k])
end
write_result(joinpath(ROOT, "results", "ritz_scan_g$tag.json"),
    (; status = "numerical, not certified", g = string(g), n0 = 6, K, Kodd,
        ed_gap_12 = gap12, ed_gap_6 = let v = sort(E6); v[2] - v[1] end,
        best = [(; k, best = best[k]) for k in ks], grid = rows,
        metadata = model_metadata(ROOT, (; g = string(g), K, Kodd, beta = "1//8:1//8:16");
            dirs = (), files = ("scripts/sectors.jl", "scripts/ritz_scan.jl"),
            packages = (CertifiedGaps,))); compact = true)
