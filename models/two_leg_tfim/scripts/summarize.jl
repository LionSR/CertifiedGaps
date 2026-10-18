# Regenerates the "## Results" section of METHOD.md from the scan, spectra (scripts/scan.jl),
# extrapolation and certificate files in results/.
#     julia --project=julia models/two_leg_tfim/scripts/summarize.jl
using CertifiedGaps.Results, Printf
const ROOT = dirname(@__DIR__)
fields = ["2_1", "5_2", "3_1", "4_1"]
result(name) = read_exact_json(joinpath(ROOT, "results", name))
gap8(tag) = (e = result("spectra_g$tag.json")["records"][end]["energies"]; e[2] - e[1])
rat(x) = Float64(parse(Rational{BigInt}, x))
const EXTRAPOLATED = Printf.Format("- `g = %s`: \$\\xi_\\delta = %.6f\$ rungs; estimated " *
    "\$(n_0, \\beta_0) = (%d, %.3f)\$, requiring a %d-rung doubled ring.\n")
const CERTIFIED = Printf.Format("- **g = %s: Δ_L ≥ %.9f** (interval bound `%s`; exact chain " *
    "bound `%s`) for every even L ≥ %s rungs. Seed (n0, β0) = (%s, %s), εₛ ≤ %s, εₜ ≤ %s; " *
    "bound / 8-rung ED gap = %.6f. Certificate: `results/certificate_g%s.json`.\n")
io = IOBuffer()
println(io, "All scan values below are **numerical, not certified**. The grid is ",
    raw"$\beta = 1/8, 2/8, \dots, 32$; $n_0 = 2, 3, 4$.", "\n")
println(io, "| g | best eligible n0 | β | numerical bound | ED gap (8 rungs) | ",
    "spectral runtime (s) |")
println(io, "|---|---:|---:|---:|---:|---:|")
for tag in fields
    c = result("scan_g$tag.json"); b = c["best"]
    t = sum(r["seconds"] for r in result("spectra_g$tag.json")["records"])
    if b === nothing
        @printf(io, "| %s | — | — | no seed | %.9f | %.2f |\n", c["g"], gap8(tag), t)
    else
        @printf(io, "| %s | %s | %.3f | %.9f | %.9f | %.2f |\n", c["g"], b["n0"], b["beta0"],
            b["bound"], gap8(tag), t)
    end
end
println(io, "\nSpectral runtimes sum the unverified block eigensolves of ",
    "`VerifiedSpectra.float_ring_spectrum` on a shared machine; they exclude Julia startup, ",
    raw"the $\beta$-grid scan, and JSON I/O.", "\n")
println(io, "### Purity and ground-energy finite-size diagnostics\n")
println(io, "| g | n | max_β min(Ps,Pτ) | δ(n)=E0(2n)−2E0(n) |")
println(io, "|---|---:|---:|---:|")
for tag in fields
    c = result("scan_g$tag.json")
    for r in c["summary"]
        @printf(io, "| %s | %s | %.9f | %.9f |\n", c["g"], r["n0"], r["max_min_purity"],
            r["delta"])
    end
end
println(io, "\nFor the fields without an ", raw"$n_0 \le 4$", " seed, a two-point fit ",
    raw"$\delta(n) = \delta(4) \exp(-(n-4)/\xi_\delta)$ to $n = 3, 4$", ", combined with a ",
    "dilute-excitation extrapolation of the 8-rung partition function, provides heuristic ",
    "estimates of the first admissible even seed (`scripts/extrapolate.jl`):\n")
for r in result("extrapolation.json")["out"]
    b = r["first_predicted_even_seed"]
    n0 = parse(Int, string(b["n0"]))
    Printf.format(io, EXTRAPOLATED, r["g"], r["xi_delta"], n0, b["beta0"], 2n0)
end
println(io, "\n### Certified results\n")
for tag in fields
    isfile(joinpath(ROOT, "results", "certificate_g$tag.json")) || continue
    c = result("certificate_g$tag.json"); s = c["seed"]; ib = c["interval_bound_decimal"]
    Printf.format(io, CERTIFIED, c["params"]["g"], floor(ib * 1e9) / 1e9, repr(ib),
        c["gap_lower_bound"], s["n0"], s["n0"], s["beta0"], rat(s["eps_s"]), rat(s["eps_t"]),
        ib / gap8(tag), tag)
end
println(io, "\nCertified: the exact RP identity, whole-spectrum eigenvalue enclosures for the ",
    "selected parameter/seed, outward log-partition and rational defect bounds, and the ",
    "resulting uniform gap under the repository's two-purity theorem. Numerical: seed ",
    "optimization, all nonverified scan spectra and finite-size gaps, the literature critical ",
    "point, and comparisons to finite-ring reference gaps. The ED reference gap is not claimed ",
    "to be an infinite-volume limit.")
p = joinpath(ROOT, "METHOD.md"); old = read(p, String)
write(p, first(split(old, "## Results")) * "## Results\n\n" * String(take!(io)))
