include("model.jl")
using Printf
root=@__DIR__
records=[JSON3.read(read(joinpath(root,"results",f),String),Dict{String,Any})
    for f in readdir(joinpath(root,"results")) if startswith(f,"summary_D")]
sort!(records;by=x->Float64(parse(Rational{Int},x["D"])))
for r in records
    d=parse(Rational{Int},r["D"])
    path=joinpath(root,"results","beta_refinement_D$(numerator(d))_$(denominator(d)).json")
    isfile(path) || continue
    fine=JSON3.read(read(path,String),Dict{String,Any})
    b=fine["best"]
    if b!==nothing && (r["best"]===nothing || b["bound"]>r["best"]["bound"])
        r["best"]=b
    end
    diag=only(x for x in r["diagnostics"] if x["n0"]==6)
    diag["max_min_purity"]=max(diag["max_min_purity"],
        maximum(min(1-x["eps_s"],1-x["eps_t"]) for x in fine["rows"]))
end
certs=[JSON3.read(read(joinpath(root,"results",f),String),Dict{String,Any})
    for f in readdir(joinpath(root,"results")) if startswith(f,"certificate_D")]
sort!(certs;by=x->Float64(parse(Rational{Int},x["D"])))
fmt(x)=@sprintf("%.9g",x)
function lower_decimal(x)
    r=Rational{BigInt}(x); n=fld(numerator(r)*big(10)^9,denominator(r))
    string(n÷big(10)^9,'.',lpad(string(n%big(10)^9),9,'0'))
end
curve=Any[]
for r in records
    d=Float64(parse(Rational{Int},r["D"])); best=r["best"]
    gap=only(x for x in r["sizes"] if x["L"]==12)["gap"]
    push!(curve,Dict("D"=>d,"bound"=>best===nothing ? nothing : best["bound"],
        "ED_gap_L12"=>gap,"ratio"=>best===nothing ? nothing : best["bound"]/gap,
        "kind"=>"Numerical scan"))
end
for c in certs
    push!(curve,Dict("D"=>Float64(parse(Rational{Int},c["D"])),"bound"=>c["bound"],
        "ED_gap_L12"=>c["ED_gap_L2n"],"ratio"=>c["ratio_to_ED_gap"],"kind"=>"Certified bound / ED gap"))
end
write_result(joinpath(root,"results","curve.json"),Dict("rows"=>curve,
    "source_hashes"=>source_hashes(root;dirs=(),files=("report.jl",)),
    "input_hashes"=>file_hashes(joinpath(root,"results"),filter(f->
        startswith(f,"summary_D") || startswith(f,"beta_refinement_D") ||
        startswith(f,"certificate_D"),readdir(joinpath(root,"results")))),
    "denominator_status"=>"Numerical finite-ring ED gap at L=12, not a certified thermodynamic gap"))
open(joinpath(root,"results","curve.csv"),"w") do io
    println(io,"D,bound,ED_gap_L12,ratio,status")
    for r in curve
        println(io,join([r["D"],something(r["bound"],"NA"),r["ED_gap_L12"],
            something(r["ratio"],"NA"),r["kind"]],','))
    end
end
spec=Dict(
    "\$schema"=>"https://vega.github.io/schema/vega/v6.json",
    "usermeta"=>Dict("source_hashes"=>source_hashes(root;dirs=(),files=("report.jl",)),
        "curve_sha256"=>sha256_file(joinpath(root,"results","curve.json"))),
    "width"=>640,"height"=>350,"padding"=>20,"background"=>"white",
    "title"=>"Spin-1 large-D: gap bound / finite-ring ED gap",
    "data"=>[Dict("name"=>"points","values"=>filter(r->r["ratio"]!==nothing,curve)),
        Dict("name"=>"scan","source"=>"points","transform"=>[
            Dict("type"=>"filter","expr"=>"datum.kind === 'Numerical scan'")])],
    "scales"=>[
        Dict("name"=>"x","type"=>"linear","domain"=>[1,4],"range"=>"width","zero"=>false),
        Dict("name"=>"y","type"=>"linear","domain"=>Dict("data"=>"points","field"=>"ratio"),"range"=>"height","nice"=>true,"zero"=>true),
        Dict("name"=>"color","type"=>"ordinal","domain"=>["Numerical scan","Certified bound / ED gap"],"range"=>["#376b9e","#b64b25"])],
    "axes"=>[Dict("orient"=>"bottom","scale"=>"x","title"=>"Single-ion anisotropy D"),
        Dict("orient"=>"left","scale"=>"y","title"=>"Bound / Δ₁₂ (numerical denominator)","grid"=>true)],
    "legends"=>[Dict("fill"=>"color","orient"=>"top-left")],
    "marks"=>[
        Dict("type"=>"line","from"=>Dict("data"=>"scan"),
            "encode"=>Dict("enter"=>Dict("x"=>Dict("scale"=>"x","field"=>"D"),
                "y"=>Dict("scale"=>"y","field"=>"ratio"),
                "stroke"=>Dict("value"=>"#376b9e"),"strokeWidth"=>Dict("value"=>2)))),
        Dict("type"=>"symbol","from"=>Dict("data"=>"points"),
            "encode"=>Dict("enter"=>Dict("x"=>Dict("scale"=>"x","field"=>"D"),
                "y"=>Dict("scale"=>"y","field"=>"ratio"),"size"=>Dict("value"=>65),
                "fill"=>Dict("scale"=>"color","field"=>"kind")))),
        Dict("type"=>"text","encode"=>Dict("enter"=>Dict("x"=>Dict("value"=>5),"y"=>Dict("value"=>330),
            "text"=>Dict("value"=>"Missing points: no successful seed on the tested grid; no gaplessness claim."),
            "fontSize"=>Dict("value"=>11),"fill"=>Dict("value"=>"#555"))))])
write_result(joinpath(root,"results","curve.vega.json"),spec)

open(joinpath(root,"REPORT.md"),"w") do io
    println(io,"# Spin-one large-D spectral-gap certificate\n")
    println(io,"Model: H = Σ Sⱼ·Sⱼ₊₁ + D Σ (Sᶻⱼ)² on periodic rings, with one site per cell.\n")
    if length(certs)<2
        println(io,"**DRAFT — full dense verification is still in progress. Only populated ",
            "certificate rows below are completed results; descriptions of the remaining ",
            "certification workflow are not claims of a completed uniform certificate.**\n")
    end
    if !isempty(certs)
        c=first(certs)
        println(io,"**Smallest certified sampled D: $(c["D"]) = $(fmt(Float64(parse(Rational{Int},c["D"])))).** ",
            "This is a finite-size, finite-grid search threshold, not a proven minimum over all real D or a phase boundary. ",
            "At that D the uniform bound is **$(lower_decimal(c["bound"]))**, for **every even L ≥ $(c["n0"]) cells = sites**.\n")
    end
    println(io,"## Certified results\n")
    println(io,"| D (exact) | n₀ | β₀ (exact) | Uniform lower bound | ED gap at 2n₀ | Bound / ED gap |")
    println(io,"|---|---:|---|---:|---:|---:|")
    for c in certs
        println(io,"| $(c["D"]) | $(c["n0"]) | $(c["beta"]) | $(lower_decimal(c["bound"])) | $(fmt(c["ED_gap_L2n"])) | $(fmt(c["ratio_to_ED_gap"])) |")
    end
    println(io,"\nDisplayed certified lower bounds are rounded down to nine decimal places; ",
        "the full binary endpoints are stored in JSON. The ED reference is the computed finite-ring L=12 gap. Its use as a thermodynamic ",
        "reference, and the plotted ratios, are numerical. Exact rational defect bounds are in ",
        "`results/certificate_D*.json`; each certificate hashes its source and full enclosed spectra.\n")
    println(io,"## Reflection positivity and rigorous method\n")
    println(io,"In the Cartesian basis, Gₓ=[0 0 0;0 0 1;0 -1 0], ",
        "Gᵧ=[0 0 -1;0 0 0;1 0 0], G_z=[0 1 0;-1 0 0;0 0 0]. ",
        "The exact identity is −h=Σₐ Gₐ⊗Gₐ, with all coefficients exactly 1, and the onsite ",
        "term D diag(1,1,0) is real symmetric. The checker independently verifies ",
        "h=P−|Ω⟩⟨Ω| and the spin commutator and Casimir.\n")
    println(io,"Every certification sector is diagonalized densely. Interval block construction ",
        "bounds matrix error; computed residuals AQ−QΛ and Q*Q−I are enlarged by rigorous ",
        "BLAS roundoff norm bounds. For residual bound r, orthogonality bound e<1, matrix ",
        "error c and norm bound M, the polar/Weyl theorem gives eigenvalue radius\n\n",
        "    η = c + r + (M + ‖Λ‖₂) e / (1 + sqrt(1−e)).\n\n",
        "[METHOD.md](METHOD.md) proves this theorem and specifies the floating-point assumptions. ",
        "Outward-rounded log-partition intervals yield exact rational upper defect bounds; ",
        "the unmodified `certified_gap_bound` then supplies the reported lower bound. ",
        "The uniform conclusion invokes the repository's stated two-purity theorem.\n")
    println(io,"## Search and numerical results\n")
    println(io,"The quarter-step grid covers D=1:0.25:4, with additional sampled D near the ",
        "first successful interval. β=0.125:0.125:20, n₀=2,…,6; selected near-threshold ",
        "samples additionally use β=2.75:0.002:3.75 at n₀=6. The table uses the best ",
        "available grid at each D. Only even n₀ enter ",
        "`best_seed`; the scope note explicitly requires even n₀. Odd n₀ are diagnostics. ",
        "The search uses Sz, translation, charge inversion, spin-inversion parity at Sz=0, ",
        "and conjugate momentum pairing. Partial Krylov spectra are permitted only for ",
        "numerical screening; their tails use the same cutoff bookkeeping as the core. ",
        "No partial spectrum is used in a certificate.\n")
    println(io,"The bound uses the core's asymmetric two-defect recurrence. It can enter ",
        "the admissible square after several iterations even when the initial spatial ",
        "defect exceeds 0.12731; that number is an equal-defect sufficient condition, ",
        "not a necessary cutoff. Thus a successful asymmetric seed can have ",
        "min(Pₛ,Pτ) below 0.87269. The certificates record the first admissible iteration.\n")
    println(io,"| D | Best n₀ | Best β₀ | Numerical bound | Δ₁₂ | Scan seconds |")
    println(io,"|---:|---:|---:|---:|---:|---:|")
    for r in records
        d=Float64(parse(Rational{Int},r["D"])); b=r["best"]
        gap=only(x for x in r["sizes"] if x["L"]==12)["gap"]
        text=b===nothing ? "— | — | no seed" : "$(b["n0"]) | $(b["beta"]) | $(fmt(b["bound"]))"
        println(io,"| $d | $text | $(fmt(gap)) | $(fmt(sum(x["seconds"] for x in r["sizes"]))) |")
    end
    println(io,"\n![Bound divided by numerical L=12 ED gap](results/curve.svg)\n")
    println(io,"At the original D=4 survey seed β=3, this numerical rerun gives ",
        "0.543863733, differing from the supplied 0.543842 by about 0.0000217. ",
        "The quoted survey value is not used as a certificate input; full dense ",
        "enclosures are computed independently. The optimized D=4 seed uses β=2.25. ",
        "The curve and its machine-readable data are in ",
        "`results/curve.svg`, `curve.csv`, and `curve.json`. Missing values mean that the ",
        "tested grid did not yield a seed, not that the physical gap vanishes.\n")
    println(io,"## Sizes, sectors, and runtimes\n")
    println(io,"Numerical rings: L=2,3,4,5,6,8,10,12. Certified spectra: L=6 and 12 at ",
        "the two stated D values. At L=12 the product dimension is 3¹²=531441. ",
        "Complete sector-by-sector dimensions, multiplicities, residual bounds and times ",
        "are stored in each `enclosed_*.json`. Each run asserts their weighted sum is 3^L.\n")
    println(io,"| D | L | Computed blocks | Largest block | Weighted dimension | Max radius | Seconds |")
    println(io,"|---|---:|---:|---:|---:|---:|---:|")
    for c in certs, L in (Int(c["n0"]),2Int(c["n0"]))
        D=parse(Rational{Int},c["D"])
        s=JSON3.read(read(joinpath(root,"results","enclosed_D$(numerator(D))_$(denominator(D))_L$L.json"),String),Dict{String,Any})
        println(io,"| $D | $L | $(length(s["sectors"])) | $(Int(maximum(x["dimension"] for x in s["sectors"]))) | $(s["dimension"]) | $(fmt(maximum(s["radii"]))) | $(fmt(s["seconds"])) |")
    end
    println(io,"\nL=14 was not attempted densely: the nonzero-charge momentum blocks have ",
        "41835 states at Sz=1, making a complex matrix alone 28.003 GB before eigenvectors ",
        "and workspaces, above this host's 16 GiB memory limit. Also n₀=7 is odd and cannot ",
        "directly invoke the stated moment theorem. Longer even seeds would require L≥16.\n")
    println(io,"## Negative-grid diagnostics and extrapolation\n")
    println(io,"All values below are numerical. For each D without a successful even seed, ",
        "the table lists max over β of min(Pₛ,Pτ), plus δ(n)=E₀(2n)−2E₀(n). ",
        "Odd n₀ diagnostics are present in the JSON summaries but do not meet the theorem.\n")
    println(io,"| D | max min purity, n=2 | n=4 | n=6 | δ(2) | δ(4) | δ(6) | Extrapolated (n,β) |")
    println(io,"|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in records
        r["best"]===nothing || continue
        d=Float64(parse(Rational{Int},r["D"])); ds=[only(x for x in r["diagnostics"] if x["n0"]==n) for n in (2,4,6)]
        gap=only(x for x in r["sizes"] if x["L"]==12)["gap"]
        ratio=ds[3]["delta"]/ds[2]["delta"]; pred="unresolved"
        if 0<ratio<1
            xi=-2/log(ratio)
            for n in 8:2:1000
                beta=log(80n)/gap
                delta=ds[3]["delta"]*exp(-(n-6)/xi)
                if beta*delta<0.075
                    pred="($n, $(fmt(beta)))"; break
                end
            end
        end
        println(io,"| $d | ",join([fmt(x["max_min_purity"]) for x in ds]," | ")," | ",
            join([fmt(x["delta"]) for x in ds]," | ")," | $pred |")
    end
    println(io,"\nExtrapolation fits δ(n)=δ(6) exp[−(n−6)/ξδ] using n=4,6 only, with ",
        "ξδ=−2/log(δ(6)/δ(4)). It holds the L=12 gap fixed and sets β=log(80n)/Δ₁₂, ",
        "then selects the first even n with βδ(n)<0.075. This corresponds to the crude ",
        "low-temperature estimates εₛ≈βδ(n)+2n exp(−βΔ₁₂)<0.1 and ",
        "ετ≈4n exp(−βΔ₁₂)=0.05. It is an illustrative size estimate, not an enclosure, ",
        "and may be poor near criticality. No extrapolated seed was certified.\n")
    println(io,"## Historical provenance and current-core replay\n")
    println(io,"The dense computations and original 123-check replay used core revision ",
        "`ef0a3e6cf54ec56dcef5e8a37963d0184f6fd31b` (`ThermalMPOGap`). Upstream then ",
        "renamed the package to `CertifiedGaps` at `22b2759e0790197b83e2571b465aa8dc69ebc0be`. ",
        "Original source files are preserved byte for byte in `provenance/legacy/`; no ",
        "historical source or core hashes were rewritten. Current executable imports were ",
        "migrated only within this certificate directory.\n")
    println(io,"The current checker verifies that the model, seed assembler and dense drivers ",
        "differ from their archived versions only by the package identifier. It verifies ",
        "historical core hashes against immutable Git objects, requires the current ",
        "recurrence implementation to be byte-identical, and replays the partition ",
        "intervals and rational bounds with the renamed core. The separate record is ",
        "`results/validation-current.json`; logs are `fast-tests.log` and ",
        "`fast-tests-legacy.log`. `check_pinned.jl` additionally replays the original ",
        "checker with the exact old core in a temporary workspace. These are replay and ",
        "migration checks, not a second full L=12 ED computation. The explicit user-requested ",
        "`certificates/spin1_largeD/` location is retained.\n")
    println(io,"## Verification and limitations\n")
    println(io,"- **Certified:** exact RP identity, complete sector counting, dense spectrum ",
        "enclosures and rational seeds for the two reported parameters, and the recurrence ",
        "lower bounds under the repository's theorem.\n",
        "- **Numerical, not certified:** all search rows, the ED reference gaps/ratios, the ",
        "curve between certified points, nonexistence of seeds at unsuccessful grid points, ",
        "and all size extrapolations.\n",
        "- No proof is claimed that all smaller real D fail, that no other β succeeds, ",
        "or that a larger even seed cannot work.\n",
        "- `test/small.jl`: independent full dense comparisons through L=6 and enclosure ",
        "checks through L=4; 20 checks passed.\n",
        "- `test/symmetry.jl`: exact transition covariance under translation and spin ",
        "inversion, charge conservation through L=6, and orbit-dimension counts through ",
        "L=12; 3278 checks passed in 4.8 seconds.\n",
        "- `test/runtests.jl`: exact RP, source/data hashes, multiplicity accounting, ",
        "partition interval reconstruction, and recurrence from stored rational seeds. ",
        "See `results/fast-tests.log` for the actual final test count and timing.\n",
        "- Julia 1.12.4, IntervalArithmetic 1.0.12, one BLAS thread. ",
        "Runtime/dependency details are in README.md and result metadata. ",
        "All project changes stay inside this directory; no core or CI changes and no merge were performed.\n")
end
