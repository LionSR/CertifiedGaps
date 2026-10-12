module LongitudinalCertificate
using LinearAlgebra, CertifiedGaps, IntervalArithmetic, JSON3
using CertifiedGaps.Seeds, CertifiedGaps.Results, CertifiedGaps.Models
BLAS.set_num_threads(1)
@assert !get_zero_subnormals()
const ROOT = @__DIR__
const GS = [1//4,1//2,3//4,1//1,3//2,2//1]
const HS = [1//20,1//10,1//4,1//2,1//1]
rat(x::AbstractString) = let a=split(x,"//"); parse(BigInt,a[1])//parse(BigInt,a[2]) end
exact(x::AbstractFloat) = Rational{BigInt}(x)
str(x::Rational) = string(numerator(x),"//",denominator(x))
str(x::AbstractFloat) = str(exact(x))

"Exact bond RP and real symmetric on-site identities, with Pauli (not spin) normalization."
function rp_check(g,h)
    X=Rational{BigInt}[0 1;1 0]; Z=Rational{BigInt}[1 0;0 -1]
    bond=Rational{BigInt}[-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]
    onsite=-g*X-h*Z
    @assert -bond == (1//1)*kron(Z,Z)
    @assert isreal(onsite) && issymmetric(onsite)
    @assert g>0 && h>=0
    true
end
rotate(s,L)=((s<<1)&((1<<L)-1))|(s>>(L-1))
"Exact integer translation orbits; offsets satisfy s=T^offset representative."
function orbits(L)
    ids=zeros(Int,1<<L); shifts=similar(ids); reps=Int[]; periods=Int[]
    for s in 0:(1<<L)-1
        ids[s+1]!=0 && continue
        push!(reps,s); t=s; p=0
        while true
            ids[t+1]=length(reps); shifts[t+1]=p; p+=1; t=rotate(t,L)
            t==s && break
        end
        push!(periods,p)
    end
    @assert sum(periods)==1<<L
    (;L,ids,shifts,reps,periods)
end
function sector_dimensions(o)
    dims=[count(p->mod(k*p,o.L)==0,o.periods) for k in 0:o.L-1]
    @assert sum(dims)==1<<o.L
    dims
end

"Construct an exactly defined Fourier block; interval mode includes basis coefficients."
function block(o,k,g,h; rigorous=false)
    L=o.L; selected=findall(p->mod(k*p,L)==0,o.periods)
    index=zeros(Int,length(o.reps)); index[selected]=eachindex(selected); d=length(selected)
    if rigorous
        iq(x)=CertifiedGaps.Recurrence.to_interval(x)
        phases=[complex(cospi(iq(2*k*t//L)),sinpi(iq(2*k*t//L))) for t in 0:L-1]
        A=fill(complex(iq(0),iq(0)),d,d)
        for (b,j) in enumerate(selected)
            s=o.reps[j]; mag=L-2count_ones(s)
            zz=L-2count_ones(xor(s,rotate(s,L)))
            A[b,b]+= -iq(zz)-iq(h)*mag
            for site in 0:L-1
                t=xor(s,1<<site); a=index[o.ids[t+1]]; a==0 && continue
                coeff=-iq(g)*sqrt(iq(o.periods[j]//o.periods[o.ids[t+1]]))
                A[a,b]+=coeff*phases[o.shifts[t+1]+1]
            end
        end
        H=zeros(ComplexF64,d,d); errors=zeros(BigFloat,d)
        setrounding(BigFloat,RoundUp) do
            for b in 1:d, a in 1:b
                z=A[a,b]; r=mid(real(z)); imv=a==b ? 0.0 : mid(imag(z))
                H[a,b]=complex(r,imv); H[b,a]=conj(H[a,b])
                er=max(BigFloat(sup(real(z)))-BigFloat(r),BigFloat(r)-BigFloat(inf(real(z))))
                ei=max(BigFloat(sup(imag(z)))-BigFloat(imv),BigFloat(imv)-BigFloat(inf(imag(z))))
                errors[a]+=er+ei
                a!=b && (errors[b]+=er+ei)
            end
        end
        return H,maximum(errors;init=big(0))
    end
    phases=[cispi(2*k*t/L) for t in 0:L-1]; H=zeros(ComplexF64,d,d)
    for (b,j) in enumerate(selected)
        s=o.reps[j]; H[b,b]=-(L-2count_ones(xor(s,rotate(s,L))))-Float64(h)*(L-2count_ones(s))
        for site in 0:L-1
            t=xor(s,1<<site); a=index[o.ids[t+1]]; a==0 && continue
            H[a,b]+=-Float64(g)*sqrt(o.periods[j]/o.periods[o.ids[t+1]])*phases[o.shifts[t+1]+1]
        end
    end
    @assert norm(H-H',Inf)<1e-10
    Matrix(Hermitian(H)),big(0)
end

"Outward Frobenius norm of a matrix whose floating entries denote exact binary numbers."
function frob_upper(A)
    setrounding(BigFloat,RoundUp) do
        s=big(0.0)
        for z in A
            a=abs(BigFloat(real(z))); b=abs(BigFloat(imag(z)))
            s+=a*a+b*b
        end
        sqrt(s)
    end
end

"Polar-factor/Weyl enclosure with explicit floating matrix-product error bounds."
function enclosed_eigen(H,assembly_error)
    d=size(H,1); @assert d<10^6; F=eigen(Hermitian(H)); Q=F.vectors; lam=F.values
    @assert maximum(abs,Q)<1e6 && maximum(abs,H)<1e6 && maximum(abs,lam)<1e6
    R=H*Q-Q.*transpose(lam); E=Q'*Q-I
    hnorm=frob_upper(H); qnorm=frob_upper(Q)
    rnorm=frob_upper(R); enorm=frob_upper(E); lnorm=BigFloat(maximum(abs,lam))
    e,r,eta=setrounding(BigFloat,RoundUp) do
        u=big(2.0)^(-53); m=16d+32; gamma=m*u/(1-m*u)
        # Factor 8 safely bounds real/imaginary dot products, subtraction, and scaling.
        tiny=big(2.0)^(-900)
        e=enorm+8gamma*(qnorm*qnorm+sqrt(BigFloat(d)))+tiny
        r=rnorm+8gamma*qnorm*(hnorm+lnorm)+tiny
        @assert e<1//2
        # ||Q-U|| <= e for e<1: avoids a downward square-root denominator.
        eta=r+(hnorm+lnorm)*e+assembly_error
        (e,r,eta)
    end
    (;values=lam,radius=2eta,orthogonality=e,residual=r,assembly_error)
end

"All momentum sectors are diagonalized densely, without Lanczos or multiplicity inference."
function ring_spectrum(L,g,h;rigorous=false)
    o=orbits(L); dims=sector_dimensions(o); E=Float64[]; records=Any[]
    start=time()
    for k in 0:L-1
        dims[k+1]==0 && continue
        H,a=block(o,k,g,h;rigorous)
        if rigorous
            enc=enclosed_eigen(H,a); append!(E,enc.values)
            push!(records,Dict("k"=>k,"dimension"=>dims[k+1],
                "values"=>str.(enc.values),"radius"=>str(enc.radius),
                "orthogonality"=>str(enc.orthogonality),"residual"=>str(enc.residual),
                "assembly_error"=>str(enc.assembly_error)))
        else
            append!(E,eigvals(Hermitian(H)))
        end
    end
    @assert length(E)==sum(dims)==1<<L
    (;E=sort!(E),dimensions=dims,blocks=records,seconds=time()-start)
end
logz(E,b)=let c=minimum(E); -b*c+log(sum(exp.(-b.*(E.-c)))) end
function defects(spec,n,b)
    a=logz(spec[n].E,b); c=logz(spec[2n].E,b); d=logz(spec[2n].E,2b)
    (-expm1(c-2a),-expm1(d-2c))
end
"Interval log partition function with a common exact shift; all eigenvalues counted."
function enclosed_logz(blocks,b)
    iq(x)=CertifiedGaps.Recurrence.to_interval(x)
    shift=minimum(rat(v) for block in blocks for v in block["values"])
    total=iq(0)
    for block in blocks
        radius=iq(rat(block["radius"]))
        for v in block["values"]
            value=iq(rat(v)-shift)+interval(-sup(radius),sup(radius))
            total+=exp(-iq(b)*value)
        end
    end
    -iq(b)*iq(shift)+log(total)
end
function seed_enclosure(small,large,b)
    a=enclosed_logz(small,b); c=enclosed_logz(large,b); d=enclosed_logz(large,2b)
    es=exact(sup(1-exp(c-2a))); et=exact(sup(1-exp(d-2c)))
    @assert 0<=es<1 && 0<=et<1
    (;es,et,logs=[a,c,d])
end
"Paths of the 30 per-point certificates (gzip-compressed JSON, `results/g*.json.gz`)."
result_paths()=sort(filter(p->startswith(basename(p),"g") && endswith(p,".json.gz"),
    readdir(joinpath(ROOT,"results");join=true)))
"Read a per-point certificate; `.gz` files are decompressed with the system `gzip`."
read_point(path)=JSON3.read(endswith(path,".gz") ? read(`gzip -dc $path`,String) :
    read(path,String),Dict{String,Any})
function hashes()
    merge(source_hashes(ROOT;dirs=("test",),files=Tuple(filter(p->endswith(p,".jl"),readdir(ROOT)))),
        Dict("../../julia/"*k=>v for (k,v) in source_hashes(joinpath(ROOT,"../../julia");
        dirs=("src",),files=("Project.toml","Manifest.toml"))))
end
function metadata(config;certified=false)
    run_metadata(config;source_sha256=hashes(),packages=(CertifiedGaps,IntervalArithmetic,JSON3),
        status=certified ? "certified" : "numerical, not certified",certified)
end
end
