include("common.jl")
include("rp.jl")
"Frobenius upper bound: IEEE sum-of-squares bound, then directed MPFR scalars."
function frob_up(A)
    count=length(A); count==0 && return BigFloat(0)
    total=sum(x->abs2(Float64(x)),A)
    @assert isfinite(total)
    setrounding(BigFloat,RoundUp) do
        k=2count+2
        gamma=BigFloat(big(k)//(big(2)^53-k))
        lower=setrounding(BigFloat,RoundDown) do; 1-gamma;end
        sqrt((BigFloat(total)+BigFloat(count)*floatmin(Float64))/lower)
    end
end
"Verified polar/Weyl enclosure radius for every ordered eigenvalue, in exact rationals."
function eigen_enclosure(A,den)
    n=size(A,1)
    H=Matrix{Float64}(A)/den
    F=eigen!(Symmetric(H)); Q=F.vectors; lam=F.values
    # Integer A is exactly represented. Verification refers to A/den, not rounded H.
    R=A*Q
    for j in 1:n,i in 1:n; R[i,j]-=Q[i,j]*lam[j]*den;end
    E=transpose(Q)*Q
    for i in 1:n; E[i,i]-=1;end
    eta,er,rr=setprecision(BigFloat,128) do
        setrounding(BigFloat,RoundUp) do
            gamma=BigFloat(big(2n+8)//(big(2)^53-(2n+8)))
            q=frob_up(Q); a=frob_up(nonzeros(A)); l=BigFloat(maximum(abs,lam))
            # Componentwise dot-product error summed in Frobenius norm.
            r=(frob_up(R)+gamma*(a*q+den*q*l)+BigFloat(n)^3*floatmin(Float64))/den
            e=frob_up(E)+gamma*(q*q+sqrt(BigFloat(n)))+BigFloat(n)^3*floatmin(Float64)
            @assert e<1//2
            sigma=setrounding(BigFloat,RoundDown) do; sqrt(1-e);end
            denom=setrounding(BigFloat,RoundDown) do; 1+sigma;end
            polar=e/denom
            radius=(r+2l*polar)/sigma
            radius,e,r
        end
    end
    radius=ceil(BigInt,eta*BigFloat(10)^12)//big(10)^12
    @assert BigFloat(radius)>=eta
    return lam,radius,Dict("dimension"=>n,"orthogonality_upper"=>string(er),"residual_upper"=>string(rr),"eigenvalue_radius"=>ratstr(radius))
end
"Recompute all Sz blocks densely; exact spin flip supplies the negative-Sz blocks."
function verified_spectrum(cells,delta)
    groups=charge_groups(cells); vals=Float64[]; radii=Rational{BigInt}[]; records=Any[]
    start=time()
    for m in 0:2cells
        A,den=integer_hamiltonian(cells,delta,groups[m+1])
        t=@elapsed lam,radius,record=eigen_enclosure(A,den)
        mult=m==0 ? 1 : 2
        record["charge"]=m;record["multiplicity"]=mult;record["seconds"]=t
        push!(records,record)
        for _ in 1:mult;append!(vals,lam);append!(radii,fill(radius,length(lam)));end
        println("VERIFIED cells=$cells M=$m dim=$(length(lam)) radius=$radius seconds=$t");flush(stdout)
    end
    @assert length(vals)==9^cells==sum(r["dimension"]*r["multiplicity"] for r in records)
    return vals,radii,records,time()-start
end
"Outward log-partition enclosure; all levels and multiplicities are included."
function verified_logz(vals,radii,beta)
    shift=minimum(vals); b=ThermalMPOGap.Recurrence.to_interval(beta)
    total=interval(0)
    for (e,r) in zip(vals,radii)
        ri=ThermalMPOGap.Recurrence.to_interval(r)
        ei=interval(e)+interval(-sup(ri),sup(ri))
        total+=exp(-b*(ei-interval(shift)))
    end
    log(total)-b*interval(shift)
end
function rational_upper(x; digits=12)
    y=Rational{BigInt}(sup(x)); d=big(10)^digits
    cld(numerator(y)*d,denominator(y))//d
end
function certify_selected(delta; n0=nothing, beta0=nothing)
    check_rp(delta)
    scan=JSON3.read(read(joinpath(ROOT,"results","scan_d$(Int(10delta)).json"),String))
    scan.best===nothing && return nothing
    n=n0===nothing ? Int(scan.best.n0) : n0;@assert n in (2,3)
    beta=beta0===nothing ? Rational{BigInt}(Float64(scan.best.beta0)) : beta0
    total_time=@elapsed begin
        v1,r1,s1,t1=verified_spectrum(n,delta)
        v2,r2,s2,t2=(2n<=4 ? verified_spectrum : verified_block_spectrum)(2n,delta)
        z1=verified_logz(v1,r1,beta);z2=verified_logz(v2,r2,beta)
        z22=verified_logz(v2,r2,2beta)
        es=rational_upper(1-exp(z2-2z1));et=rational_upper(1-exp(z22-2z2))
        @assert 0<es<1 && 0<et<1
        bound=certified_gap_bound(es,et,beta)
        @assert bound.bound!==nothing && bound.bound>0
    end
    interval_record(x)=Dict("lower"=>ratstr(Rational{BigInt}(inf(x))),"upper"=>ratstr(Rational{BigInt}(sup(x))))
    data=Dict("metadata"=>provenance(Dict("delta"=>ratstr(delta),"n0"=>n,"beta0"=>ratstr(beta));certified=true),"delta"=>ratstr(delta),"n0"=>n,"beta0"=>ratstr(beta),"epsilon_s"=>ratstr(es),"epsilon_t"=>ratstr(et),"bound"=>bound.bound,"bound_exact_binary"=>ratstr(Rational{BigInt}(bound.bound)),"logZ_n_beta"=>interval_record(z1),"logZ_2n_beta"=>interval_record(z2),"logZ_2n_2beta"=>interval_record(z22),"sectors_n"=>s1,"sectors_2n"=>s2,"seconds"=>total_time,"reference_gap_numerical"=>sort(v2)[2]-minimum(v2),"validity"=>"even L >= $n cells; site counts N divisible by 4 and N >= $(4cld(n,2))","positive_transfer_gram_proof"=>true,"verification"=>"dense symmetry sectors; exact projectors for large rings; polar/Weyl and rigorous IEEE norm bounds")
    write_result(joinpath(ROOT,"results","certificate_d$(Int(10delta)).json"),data)
    println("CERTIFICATE delta=$delta bound=$(bound.bound) seed=($n,$beta) es=$es et=$et")
    data
end
include("verified_blocks.jl")
if abspath(PROGRAM_FILE)==@__FILE__
    certify_selected(parse(Int,ARGS[1])//10; n0=length(ARGS)>1 ? parse(Int,ARGS[2]) : nothing,beta0=length(ARGS)>2 ? parserat(ARGS[3]) : nothing)
end
