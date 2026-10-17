# Regenerates the "## Results" section of METHOD.md from the scan, spectra (scripts/scan.jl),
# extrapolation and certificate files in results/.
using CertifiedGaps.Results, Printf
root=dirname(@__DIR__); fields=["2_1","5_2","3_1","4_1"]
io=IOBuffer()
println(io,"All scan values below are **numerical, not certified**. The grid is β=1/8,2/8,…,32; n0=2,3,4.\n")
println(io,"| g | best eligible n0 | β | numerical bound | ED gap (8 rungs) | spectral runtime (s) |")
println(io,"|---|---:|---:|---:|---:|---:|")
for tag in fields
    c=read_exact_json(joinpath(root,"results","scan_g$tag.json"))
    s=read_exact_json(joinpath(root,"results","spectra_g$tag.json"))
    b=c["best"]; e=s["records"][end]["energies"]
    t=sum(r["seconds"] for r in s["records"])
    if b===nothing
        @printf(io,"| %s | — | — | no seed | %.9f | %.2f |\n",c["g"],e[2]-e[1],t)
    else
        @printf(io,"| %s | %s | %.3f | %.9f | %.9f | %.2f |\n",c["g"],b["n₀"],b["β₀"],b["bound"],e[2]-e[1],t)
    end
end
println(io,"\nSpectral runtimes sum the dense-spectrum calls of the core block reduction; they exclude Julia startup, the β-grid scan, and JSON I/O.\n")
println(io,"\n### Purity and ground-energy finite-size diagnostics\n")
println(io,"| g | n | max_β min(Ps,Pτ) | δ(n)=E0(2n)−2E0(n) |")
println(io,"|---|---:|---:|---:|")
for tag in fields
    c=read_exact_json(joinpath(root,"results","scan_g$tag.json"))
    for r in c["summary"]
        @printf(io,"| %s | %s | %.9f | %.9f |\n",c["g"],r["n"],r["max_min_purity"],r["delta"])
    end
end
println(io,"\nFor the unsuccessful fields, a two-point fit δ(n)=δ(4)exp(−(n−4)/ξδ) to n=3,4, combined with an extensive dilute-excitation extrapolation of the eight-rung partition function, predicts the following first even seeds. These are **uncontrolled estimates**, not computed larger-ring spectra or certificates; near-critical finite-size effects can bias them substantially. `scripts/extrapolate.jl` records the explicit ansatz and grid.\n")
ex=read_exact_json(joinpath(root,"results","extrapolation.json"))
for r in ex["out"]
    b=r["first_predicted_even_seed"]
    @printf(io,"- g=%s: ξδ=%.6f rungs; predicted (n0,β)=(%s,%.3f), requiring %d rungs for direct verification.\n",r["g"],r["xi_delta"],b["n₀"],b["β₀"],2parse(Int,b["n₀"]))
end
println(io,"\n### Certified results\n")
for tag in fields
    p=joinpath(root,"results","certificate_g$tag.json"); isfile(p)||continue
    c=read_exact_json(p); s=c["seed"]
    e=read_exact_json(joinpath(root,"results","spectra_g$tag.json"))["records"][end]["energies"]
    ib=c["interval_bound_decimal"]
    @printf(io,"- **g=%s: Δ_L ≥ %.9f** (interval bound `%s`; exact chain bound `%s`) for every %s. Seed (n0,β)=(%s,%s), εₛ ≤ %s, εₜ ≤ %s; bound / 8-rung ED gap = %.6f. Certificate: `results/certificate_g%s.json`.\n",
        c["point"]["g"],floor(ib*1e9)/1e9,repr(ib),c["gap_lower_bound"],c["point"]["lengths"],
        s["n0"],s["beta0"],Float64(parse(Rational{BigInt},s["eps_s"])),
        Float64(parse(Rational{BigInt},s["eps_t"])),ib/(e[2]-e[1]),tag)
end
println(io,"\nCertified: the exact RP identity, whole-spectrum eigenvalue enclosures for the selected parameter/seed, outward log-partition and rational defect bounds, and the resulting uniform gap under the repository's two-purity theorem. Numerical: seed optimization, all nonverified scan spectra and finite-size gaps, the literature critical point, and comparisons to finite-ring reference gaps. The ED reference gap is not claimed to be an infinite-volume limit.\n")
p=joinpath(root,"METHOD.md"); old=read(p,String)
write(p,first(split(old,"## Results"))*"## Results\n\n"*String(take!(io)))
