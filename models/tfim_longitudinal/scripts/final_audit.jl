# Independent checks of stored artifacts: no calls to seed_enclosure or orbits.
# Usage: julia --project=julia models/tfim_longitudinal/scripts/final_audit.jl (~30 s).
using LinearAlgebra, IntervalArithmetic, CertifiedGaps, CertifiedGaps.Results, JSON3
using Test
BLAS.set_num_threads(1)
const AUDIT_ROOT = normpath(joinpath(@__DIR__, ".."))
parseq(s) = let a=split(s,"//"); parse(BigInt,a[1])//parse(BigInt,a[2]) end
qstr(q) = string(numerator(q),"//",denominator(q))
biginterval(q::Rational) = interval(BigFloat,numerator(q))/interval(BigFloat,denominator(q))

function mobius(n)
    result=1; p=2
    while p*p<=n
        if n%p==0
            n÷=p; result=-result
            n%p==0 && return 0
        end
        p+=1
    end
    n>1 ? -result : result
end
ramanujan(q,k) = sum(d*mobius(q÷d) for d in 1:gcd(q,k) if gcd(q,k)%d==0)
function character_dimensions(L)
    [begin
        numerator=sum(big(2)^(L÷q)*ramanujan(q,k) for q in 1:L if L%q==0)
        numerator%L==0 || error("character sum not divisible by L=$L")
        Int(numerator÷L)
    end for k in 0:L-1]
end

function independent_logz(blocks,beta)
    center=minimum(parseq(v) for block in blocks for v in block["values"])
    total=interval(BigFloat,0)
    for block in blocks
        radius=parseq(block["radius"])
        for v in block["values"]
            q=parseq(v)-center
            lo=inf(biginterval(q-radius)); hi=sup(biginterval(q+radius))
            total+=exp(-biginterval(beta)*interval(BigFloat,lo,hi))
        end
    end
    -biginterval(beta)*biginterval(center)+log(total)
end

started=time(); rows=Any[]; pairs=Set{Tuple{Rational{BigInt},Rational{BigInt}}}()
paths=sort(filter(p->startswith(basename(p),"certificate_g") && endswith(p,".json.gz"),
    readdir(joinpath(AUDIT_ROOT,"results");join=true)))
@testset "Independent final artifact audit" begin
    @test length(paths)==30
    for path in paths
        data=JSON3.read(read(`gzip -dc $path`,String),Dict{String,Any})
        g=parseq(data["g"]); h=parseq(data["h"]); push!(pairs,(g,h))
        cert=data["certificate"]; n=cert["n"]; beta=parseq(cert["beta"])
        @test n==6 && iseven(n) && beta>0 && g>0 && h>0
        # Independently specified Pauli matrices and exact local terms.
        Z=diagm(Rational{BigInt}[1,-1]); X=Rational{BigInt}[0 1;1 0]
        bond=diagm(Rational{BigInt}[-1,1,1,-1]); onsite=Rational{BigInt}[-h -g;-g h]
        @test -bond==kron(Z,Z) && onsite==-g*X-h*Z && issymmetric(onsite)
        for (name,L) in (("small",n),("large",2n))
            blocks=cert[name]; dims=character_dimensions(L)
            @test [b["k"] for b in blocks]==collect(0:L-1)
            @test [b["dimension"] for b in blocks]==dims
            @test sum(dims)==2^L
            tr=interval(BigFloat,0); tr2=interval(BigFloat,0)
            for block in blocks
                @test block["dimension"]==length(block["values"])
                radius=parseq(block["radius"])
                @test radius>0 && 0<=parseq(block["orthogonality"])<1//2
                for token in block["values"]
                    v=parseq(token)
                    x=interval(BigFloat,inf(biginterval(v-radius)),
                        sup(biginterval(v+radius)))
                    tr+=x; tr2+=x*x
                end
            end
            # Distinct nonidentity Pauli words are Hilbert-Schmidt orthogonal.
            expected2=big(2)^L*(L==2 ? 4+2g*g+2h*h : L*(1+g*g+h*h))
            @test inf(tr)<=0<=sup(tr)
            @test inf(tr2)<=expected2<=sup(tr2)
        end
        a=independent_logz(cert["small"],beta)
        b=independent_logz(cert["large"],beta)
        c=independent_logz(cert["large"],2beta)
        # Endpoint directions are explicit and independent of the production routine.
        lower_ps=exp(interval(BigFloat,inf(b))-2interval(BigFloat,sup(a)))
        lower_pt=exp(interval(BigFloat,inf(c))-2interval(BigFloat,sup(b)))
        es=interval(BigFloat,1)-lower_ps; et=interval(BigFloat,1)-lower_pt
        stored_s=parseq(cert["epsilon_s"]); stored_t=parseq(cert["epsilon_t"])
        @test Rational{BigInt}(sup(es))<=stored_s<1
        @test Rational{BigInt}(sup(et))<=stored_t<1
        answer=certified_gap_bound(stored_s,stored_t,beta)
        @test answer.bound!==nothing && answer.bound>0
        @test Rational{BigInt}(answer.bound)==parseq(cert["bound_exact_binary"])
        @test answer.j_admissible!==nothing
        push!(rows,Dict("g"=>qstr(g),"h"=>qstr(h),"n"=>n,"beta"=>qstr(beta),
            "independent_epsilon_s_upper"=>qstr(Rational{BigInt}(sup(es))),
            "independent_epsilon_t_upper"=>qstr(Rational{BigInt}(sup(et))),
            "stored_bound_exact"=>cert["bound_exact_binary"],
            "j_admissible"=>answer.j_admissible,"certificate_sha256"=>sha256_file(path)))
        println("Audited ",basename(path)); flush(stdout)
    end
    expected=Set((Rational{BigInt}(g),Rational{BigInt}(h))
        for g in (1//4,1//2,3//4,1//1,3//2,2//1)
        for h in (1//20,1//10,1//4,1//2,1//1))
    @test pairs==expected
end
worst=argmin(r->parseq(r["stored_bound_exact"]),rows)
write_result(joinpath(AUDIT_ROOT,"results","final_audit.json"),Dict(
    "status"=>"passed","parameters"=>rows,"worst"=>worst,"seconds"=>time()-started,
    "method"=>"Independent character dimensions, Pauli trace moments, " *
        "256-bit interval partition sums",
    "source_hashes"=>source_hashes(AUDIT_ROOT;dirs=("src","scripts","test"),
        files=Tuple(filter(p->endswith(p,".jl"),readdir(AUDIT_ROOT)))),
    "julia"=>string(VERSION),"bigfloat_precision"=>precision(BigFloat)))
