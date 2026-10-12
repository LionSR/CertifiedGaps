include("ladder.jl")
using JSON3

"Enclose log Z by interval exponentials of every eigenvalue enclosure."
function enclosed_logZ(s,beta)
    shift=s.energies[1]
    z=interval(0.0)
    for (e,r) in zip(s.energies,s.radii)
        z+=exp(-to_interval(beta)*((interval(e)+interval(-r,r))-shift))
    end
    log(z)-to_interval(beta)*shift
end

function certify()
    start=time(); J=4//1; @assert exact_rp(J)
    scanpath=joinpath(@__DIR__,"results","scan_J4p0.json")
    if isfile(scanpath)
        best=JSON3.read(read(scanpath,String)).best
    else
        specs=Dict{Int,PartialSpectrum}()
        for L in (2,4,8)
            path=joinpath(@__DIR__,"results","spectrum_J4p0_L$(L).json")
            data=JSON3.read(read(path,String))
            specs[L]=PartialSpectrum(Float64.(data.energies),4^L)
        end
        candidate=best_seed(scan_seeds(spectral_defects(specs),(2,4),0.25:0.25:16))
        best=candidate===nothing ? nothing : (;n0=candidate.n₀,beta=candidate.β₀)
    end
    @assert best !== nothing
    n=Int(best.n0); beta=Rational{BigInt}(Float64(best.beta))
    @assert iseven(n)
    a=full_spectrum(n,J;verified=true); b=full_spectrum(2n,J;verified=true)
    la=enclosed_logZ(a,beta); lb=enclosed_logZ(b,beta); lc=enclosed_logZ(b,2beta)
    es=-expm1(lb-2la); et=-expm1(lc-2lb)
    # Round rational upper bounds to a convenient decimal grid, still outward.
    ratup(x)=ceil(BigInt,Rational{BigInt}(sup(x))*big(10)^12)//big(10)^12
    eps_s=ratup(es); eps_t=ratup(et)
    answer=certified_gap_bound(eps_s,eps_t,beta); @assert answer.bound !== nothing
    enc(x)=Dict("lo"=>string(Rational{BigInt}(inf(x))),"hi"=>string(Rational{BigInt}(sup(x))))
    for (L,s) in ((n,a),(2n,b))
        write_result(joinpath(@__DIR__,"results","verified_L$(L).json"),Dict(
            "metadata"=>provenance(Dict("J"=>string(J),"L"=>L);certified=true),
            "L"=>L,"dimension"=>s.dimension,"seconds"=>s.seconds,"sectors"=>s.rows,
            "energies"=>string.(Rational{BigInt}.(s.energies)),
            "radii"=>string.(Rational{BigInt}.(s.radii))))
    end
    gap=b.energies[2]-b.energies[1]
    write_result(joinpath(@__DIR__,"results","certificate.json"),Dict(
        "metadata"=>provenance(Dict("J"=>string(J));certified=true),
        "J"=>string(J),"n0"=>n,"beta"=>string(beta),"eps_s"=>string(eps_s),
        "eps_t"=>string(eps_t),"bound"=>answer.bound,"recurrence"=>answer,
        "logZ_n_beta"=>enc(la),"logZ_2n_beta"=>enc(lb),"logZ_2n_2beta"=>enc(lc),
        "validity"=>"every even L >= $n rungs (N=2L >= $(2n) spin sites; N divisible by 4)",
        "seconds"=>time()-start,"ed_gap_2n"=>gap,"bound_over_ed_gap"=>answer.bound/gap,
        "maximum_eigenvalue_radius"=>maximum(b.radii)))
    println("CERTIFIED ",answer," eps=",eps_s," ",eps_t," beta=",beta)
end
abspath(PROGRAM_FILE)==(@__FILE__) && certify()
