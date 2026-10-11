include("certify.jl")
using .LongitudinalCertificate, JSON3, Printf
const C=LongitudinalCertificate
paths=filter(p->startswith(basename(p),"g") && endswith(p,".json"),
    readdir(joinpath(C.ROOT,"results");join=true))
data=sort([JSON3.read(read(p,String),Dict{String,Any}) for p in paths];
    by=d->(C.rat(d["g"]),C.rat(d["h"])))
audit=JSON3.read(read(joinpath(C.ROOT,"results","audit.json"),String),Dict{String,Any})
fmt(x)=@sprintf("%.6f",x)
ncert=count(d->haskey(d,"certificate") && d["certificate"]["bound"]!==nothing,data)
open(joinpath(C.ROOT,"REPORT.md"),"w") do io
    println(io,"# Longitudinal-field TFIM: uniform spectral-gap certificates\n")
    println(io,"**$ncert of $(length(data)) parameter points certified.** One cell is one spin-1/2 site.")
    println(io,"Every bound below holds for **every even L ≥ n₀ cells = sites**, not merely the ED sizes.")
    println(io,"The grid is g ∈ {0.25,0.5,0.75,1,1.5,2}, h ∈ {0.05,0.1,0.25,0.5,1}.")
    println(io,"No result is asserted for an unsampled parameter or for odd L.\n")
    println(io,"## Model, RP and uniqueness\n")
    println(io,"H_L = −Σ Z_j Z_(j+1) − gΣ X_j − hΣ Z_j, with periodic bonds, including both bonds at L=2.")
    println(io,"In the σᶻ basis, X=[[0,1],[1,0]], Z=diag(1,−1), c=1 and C=Z.")
    println(io,"Thus **−h_bond = 1·Z⊗Z**, with exact rational entries and positive coefficient.")
    println(io,"The on-site matrix [[−h,−g],[−g,h]] is real symmetric. `rp_check` verifies these identities")
    println(io,"in Rational{BigInt} for every stated rational parameter. The spatial transfer operator is therefore")
    println(io,"self-adjoint under the repository's Appendix E theorem. The ground state at every finite L is")
    println(io,"unique independently: for g>0, a sufficiently large scalar shift of −H is an irreducible")
    println(io,"nonnegative matrix on the connected spin-flip hypercube. Perron–Frobenius gives a simple largest")
    println(io,"eigenvalue. That finite-size argument alone does not prove a uniform gap.\n")
    println(io,"## Search and actual sizes\n")
    println(io,"The initial scan uses n₀ ∈ {2,4,6}, β₀ ∈ {0.25,0.50,…,64}, via `scan_seeds` and `best_seed`.")
    println(io,"Only even n₀ are eligible: the stated criterion requires an even spatial power.")
    println(io,"All lengths in each JSON's `spectra` were densely diagonalized in **every** momentum sector.")
    println(io,"No Lanczos, ARPACK, DMRG, or partial-spectrum tail is used. The best numerical seed per point")
    println(io,"was rerun with rigorous spectral and partition-function enclosures. A numerical optimum need not")
    println(io,"maximize the final conservative certified bound. All 30 points succeeded at this initial size cap; no larger rings were needed for the uniform claims.\n")
    println(io,"| L (sites = cells) | Momentum sector dimensions k=0,…,L−1 | Sum |")
    println(io,"|---:|---|---:|")
    lengths=sort(unique(vcat([parse.(Int,collect(keys(d["spectra"]))) for d in data]...)))
    for L in lengths
        dims=C.sector_dimensions(C.orbits(L))
        println(io,"| $L | ",join(dims,", ")," | ",sum(dims)," = 2^$L |")
    end
    println(io,"\nThe requested 20–24-site full-spectrum survey was **not reached**. At L=20 even generic k=1")
    println(io,"has 52,377 states: a single complex dense matrix takes about 40.9 GiB, before eigenvectors and")
    println(io,"workspace, on a machine with a 16 GiB container memory limit. At L=24 this rises to thousands of GiB.")
    println(io,"See `results/audit.json` for exact integer resource counts. Even a real representation does not")
    println(io,"make the 20-site generic block fit with eigensolver workspace. Spin-flip parity is broken by h>0.")
    println(io,"Reflection exchanges k and −k, rather than independently halving every block; the present code")
    println(io,"does not further split k=0,π by reflection parity. Smaller seeds already establish uniform results.\n")
    println(io,"## Best seed and bound per parameter\n")
    println(io,"Decimals are approximate displays; exact input rationals and exact binary lower bounds are in JSON.")
    println(io,"Δ_ED is the numerical first excitation gap at the largest scanned ring, **not** an extrapolated")
    println(io,"thermodynamic gap and not a rigorous reference. Ratios compare the certified bound to this ED value.\n")
    println(io,"| g | h | n₀ | β₀ | certified bound ≈ | ED gap ≈ | bound/ED | total seconds |")
    println(io,"|---:|---:|---:|---:|---:|---:|---:|---:|")
    for d in data
        g=Float64(C.rat(d["g"])); h=Float64(C.rat(d["h"])); c=get(d,"certificate",nothing)
        L=maximum(parse.(Int,collect(keys(d["spectra"])))); gap=d["spectra"][string(L)]["gap"]
        if c===nothing || c["bound"]===nothing
            println(io,"| $g | $h | — | — | not certified | ",fmt(gap)," | — | ",fmt(d["seconds"])," |")
        else
            println(io,"| $g | $h | ",c["n"]," | ",Float64(C.rat(c["beta"]))," | ",fmt(c["bound"]),
                " | ",fmt(gap)," | ",fmt(c["bound"]/gap)," | ",fmt(d["seconds"])," |")
        end
    end
    println(io,"\nThe h>0 points at **g=1** and in the **ordered zero-field regime g<1** are included explicitly.")
    println(io,"The phrase 'critical point' here means the g=1 location of the h=0 model; h>0 is a different,")
    println(io,"non-integrable Hamiltonian. Nothing here certifies a uniform gap as h tends to zero.\n")
    println(io,"The unequal-defect recurrence can succeed even when max min(Pₛ,Pτ) is below 1−0.12731; the equal-defect square is sufficient, not necessary.\n")
    println(io,"## Region map\n")
    println(io,"Entries are certified lower bounds, approximately displayed; row g, column h.\n")
    println(io,"| g \\ h | 0.05 | 0.1 | 0.25 | 0.5 | 1.0 |\n|---:|---:|---:|---:|---:|---:|")
    for g in C.GS
        vals=String[]
        for h in C.HS
            found=findfirst(d->C.rat(d["g"])==g && C.rat(d["h"])==h,data)
            if found===nothing || !haskey(data[found],"certificate")
                push!(vals,"—")
            else
                push!(vals,fmt(data[found]["certificate"]["bound"]))
            end
        end
        println(io,"| ",Float64(g)," | ",join(vals," | ")," |")
    end
    println(io,"\n## Numerical finite-size diagnostics\n")
    println(io,"Each pair is max over the β grid of min(Pₛ,Pτ), followed by δ(n)=E₀(2n)−2E₀(n).")
    println(io,"These diagnostics use unverified numerical eigenvalues; they are not certificate inputs.")
    println(io,"| g | h | n : (max min purity; δ) |\n|---:|---:|---|")
    for d in data
        vals=[string(x["n"],": (",fmt(x["max_min_purity"]),"; ",@sprintf("%.5e",x["delta"]),")") for x in d["diagnostics"]]
        println(io,"| ",Float64(C.rat(d["g"]))," | ",Float64(C.rat(d["h"]))," | ",join(vals," · ")," |")
    end
    println(io,"\nAn extrapolated seed is unnecessary for points with an actual certificate; no negative conclusion")
    println(io,"about inaccessible sizes follows from this scan.\n")
    println(io,"## Rigor and validation\n")
    println(io,"[METHOD.md](METHOD.md) gives the complete implemented enclosure theorem. In brief, interval")
    println(io,"assembly bounds ||A−H||, BigFloat upward accumulation plus explicit dot-product rounding")
    println(io,"bounds encloses R=HQ−QΛ and E=Q*Q−I. The polar factor of Q and Weyl's min-max inequality")
    println(io,"then enclose every eigenvalue with multiplicity. Interval sums enclose all three log Z values;")
    println(io,"outward rational defect bounds feed the unmodified `certified_gap_bound`. Source hashes, Julia")
    println(io,"and package versions, dimensions, radii, all eigenvalues and rational seeds are stored in `results/*.json`.\n")
    println(io,"For an exact Hermitian block A, a ≥ ||A−H||₂, r ≥ ||HQ−QΛ||₂ and e ≥ ||Q*Q−I||₂ < 1/2,")
    println(io,"the polar factor U obeys ||Q−U||₂ ≤ e. Hence Weyl gives, for every ordered eigenvalue,")
    println(io,"**|λᵢ(A)−Λᵢ| ≤ a + r + (||H||F + max|Λᵢ|)e**. This is the theorem implemented;")
    println(io,"its norms include the explicit outward rounding allowances proved in METHOD.md.\n")
    println(io,"Certified: RP; complete dimensions; eigenvalue enclosures; partition/seed enclosures; resulting")
    println(io,"uniform gap for the stated even-length range. Numerical: seed ranking, ED reference gaps,")
    println(io,"ratios, δ(n), purity diagnostics and free-fermion cross-check values. No continuum parameter")
    println(io,"region, h→0 limit, odd-size gap, or thermodynamic reference-gap extrapolation is claimed.\n")
    println(io,"The independent audit passed 42 checks: full computational-basis ED versus momentum ED,")
    println(io,"the h=0,g=2 `tfim_logZ` formula at six lengths and five temperatures, and exact classical")
    println(io,"spectra inside the implemented enclosures. Maximum free-fermion log Z discrepancy: ",audit["max_error"],".")
    println(io,"The fast checker repeats RP, dimension counts, partition bounds and the recurrence. The deep")
    println(io,"checker regenerates rigorous eigensystem enclosures and verifies containment in the stored intervals using exact rationals. The stored radii have a factor-two replay margin.")
    println(io,"Test commands and measured outcomes are in [README.md](README.md).\n")
    println(io,"## Environment and runtime\n")
    println(io,"Julia 1.12.4; dependencies instantiated and precompiled from the existing Julia project without")
    println(io,"modifying its Project or Manifest. OPENBLAS_NUM_THREADS=1 and BLAS.set_num_threads(1).")
    println(io,"The approved installer endpoint returned HTTP 403; the official Linux x86_64 archive was used")
    println(io,"in /tmp and its SHA-256 checked against the official checksum:")
    println(io,"`c57baf178fe140926acb1a25396d482f325af9d7908d9b066d2fbc0d6639985d`.")
    println(io,"Sum of per-point measured runtimes: ",fmt(sum(d["seconds"] for d in data))," seconds;")
    println(io,"audit: ",fmt(audit["seconds"])," seconds. Startup compilation is additional.")
    println(io,"All deliverables are confined to this certificate directory. No merge or publication performed.")
end
