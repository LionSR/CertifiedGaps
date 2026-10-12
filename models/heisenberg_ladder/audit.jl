# Independent arithmetic checks after actual full spectral replay; saved radii are not used.
include("certify.jl")
using Test

function big_logZ(s,beta)
    setprecision(BigFloat,256) do
        B=BigFloat(beta); @assert Rational{BigInt}(B)==beta
        lo=BigFloat(0); hi=BigFloat(0)
        for (e,r) in zip(s.energies,s.radii)
            lower=setrounding(BigFloat,RoundDown) do
                BigFloat(e)-BigFloat(r)
            end
            upper=setrounding(BigFloat,RoundUp) do
                BigFloat(e)+BigFloat(r)
            end
            lo=setrounding(BigFloat,RoundDown) do
                lo+exp(-B*upper)
            end
            hi=setrounding(BigFloat,RoundUp) do
                hi+exp(-B*lower)
            end
        end
        a=setrounding(BigFloat,RoundDown) do; log(lo); end
        b=setrounding(BigFloat,RoundUp) do; log(hi); end
        a,b
    end
end

function exact_trace_checks(s,L,J)
    low=big(0)//1; high=big(0)//1; low2=big(0)//1; high2=big(0)//1
    for (e,r) in zip(s.energies,s.radii)
        a=Rational{BigInt}(e)-Rational{BigInt}(r)
        b=Rational{BigInt}(e)+Rational{BigInt}(r)
        low+=a; high+=b
        low2+=a<=0<=b ? 0 : min(a^2,b^2)
        high2+=max(a^2,b^2)
    end
    expected=big(4)^L * (3//16) * L * (J^2+2)
    @test low<=0<=high
    @test low2<=expected<=high2
    Dict("L"=>L,"dimension"=>s.dimension,"trace_lower"=>string(low),
        "trace_upper"=>string(high),"trace2_lower"=>string(low2),
        "trace2_upper"=>string(high2),"trace2_exact"=>string(expected))
end

function audit()
    start=time(); c=JSON3.read(read(joinpath(@__DIR__,"results/certificate.json"),String))
    q(s)=parse(Rational{BigInt},s)
    n=Int(c.n0); beta=q(c.beta); J=q(c.J)
    @test iseven(n) && n==4 && beta==2 && J==4
    @test exact_rp(J)
    a=full_spectrum(n,4//1;verified=true)
    b=full_spectrum(2n,4//1;verified=true)
    traces=[exact_trace_checks(a,n,J),exact_trace_checks(b,2n,J)]
    la=big_logZ(a,beta); lb=big_logZ(b,beta); lc=big_logZ(b,2beta)
    ps=setrounding(BigFloat,RoundDown) do; exp(lb[1]-2la[2]); end
    pt=setrounding(BigFloat,RoundDown) do; exp(lc[1]-2lb[2]); end
    es=setrounding(BigFloat,RoundUp) do; 1-ps; end
    et=setrounding(BigFloat,RoundUp) do; 1-pt; end
    @test Rational{BigInt}(es)<=q(c.eps_s)
    @test Rational{BigInt}(et)<=q(c.eps_t)
    @test certified_gap_bound(q(c.eps_s),q(c.eps_t),beta).bound==c.bound
    enc(x)=Dict("lo"=>string(Rational{BigInt}(x[1])),"hi"=>string(Rational{BigInt}(x[2])))
    write_result(joinpath(@__DIR__,"results/audit.json"),Dict(
        "metadata"=>provenance(Dict("arithmetic"=>"256-bit directed MPFR and exact rational trace moments");certified=true),
        "stored_radii_used"=>false,"fresh_dimensions"=>[a.dimension,b.dimension],
        "traces"=>traces,"logZ_n_beta"=>enc(la),"logZ_2n_beta"=>enc(lb),
        "logZ_2n_2beta"=>enc(lc),"eps_s_upper"=>string(Rational{BigInt}(es)),
        "eps_t_upper"=>string(Rational{BigInt}(et)),"seconds"=>time()-start,
        "limits"=>"Fresh eigensolves share the symmetry builder and norm-bound implementation; MPFR partition arithmetic and exact trace moment checks are separate. The repository theorem and IEEE/OpenBLAS arithmetic model remain assumptions."))
    println("AUDIT PASSED: fresh spectra, exact trace moments, independent directed MPFR purities.")
end
audit()
