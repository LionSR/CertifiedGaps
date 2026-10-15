"""
    LongitudinalCertificate

Seed certificates for the longitudinal-field TFIM ring
`H_L = -Σ Z_j Z_{j+1} - g Σ X_j - h Σ Z_j` (Pauli matrices, one site per cell):
exact momentum-sector assembly, polar/Weyl eigenvalue enclosures with explicit rounding
allowances (see `METHOD.md`), interval partition sums and result I/O.
"""
module LongitudinalCertificate
using LinearAlgebra, CertifiedGaps, IntervalArithmetic, JSON3
using CertifiedGaps.Seeds, CertifiedGaps.Results, CertifiedGaps.Models
BLAS.set_num_threads(1)
# The rounding allowances assume gradual underflow (see METHOD.md).
get_zero_subnormals() && error("flush-to-zero is enabled; the enclosures assume gradual underflow")
"The model directory (`models/tfim_longitudinal`); results live in `ROOT/results`."
const ROOT = normpath(joinpath(@__DIR__, ".."))
"Transverse fields g of the certified grid."
const GS = [1//4,1//2,3//4,1//1,3//2,2//1]
"Longitudinal fields h of the certified grid."
const HS = [1//20,1//10,1//4,1//2,1//1]
"Parse an exact `\"num//den\"` string."
rat(x::AbstractString) = let a=split(x,"//"); parse(BigInt,a[1])//parse(BigInt,a[2]) end
"The exact rational value of a binary floating-point number."
exact(x::AbstractFloat) = Rational{BigInt}(x)
"Format an exact number as a `\"num//den\"` string."
str(x::Rational) = string(numerator(x),"//",denominator(x))
str(x::AbstractFloat) = str(exact(x))

"""
    @check cond

Throw an error naming `cond` unless it holds. Used instead of `@assert`, which Julia may
elide; every condition here is part of the proof and must always run.
"""
macro check(ex)
    msg = string(ex)
    return :($(esc(ex)) || error("check failed: " * $msg))
end

"Exact bond RP and real symmetric on-site identities, with Pauli (not spin) normalization."
function rp_check(g,h)
    X=Rational{BigInt}[0 1;1 0]; Z=Rational{BigInt}[1 0;0 -1]
    bond=Rational{BigInt}[-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]
    onsite=-g*X-h*Z
    @check -bond == (1//1)*kron(Z,Z)
    @check isreal(onsite) && issymmetric(onsite)
    @check g>0 && h>=0
    true
end
"Cyclic translation of the L-bit configuration `s` by one site."
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
    @check sum(periods)==1<<L
    (;L,ids,shifts,reps,periods)
end
"Dimensions of the momentum sectors k = 0, …, L-1 for the orbits `o`."
function sector_dimensions(o)
    dims=[count(p->mod(k*p,o.L)==0,o.periods) for k in 0:o.L-1]
    @check sum(dims)==1<<o.L
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
    @check norm(H-H',Inf)<1e-10
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
    d=size(H,1); @check d<10^6; F=eigen(Hermitian(H)); Q=F.vectors; lam=F.values
    @check maximum(abs,Q)<1e6 && maximum(abs,H)<1e6 && maximum(abs,lam)<1e6
    R=H*Q-Q.*transpose(lam); E=Q'*Q-I
    hnorm=frob_upper(H); qnorm=frob_upper(Q)
    rnorm=frob_upper(R); enorm=frob_upper(E); lnorm=BigFloat(maximum(abs,lam))
    e,r,eta=setrounding(BigFloat,RoundUp) do
        u=big(2.0)^(-53); m=16d+32; gamma=m*u/(1-m*u)
        # Factor 8 safely bounds real/imaginary dot products, subtraction, and scaling.
        tiny=big(2.0)^(-900)
        e=enorm+8gamma*(qnorm*qnorm+sqrt(BigFloat(d)))+tiny
        r=rnorm+8gamma*qnorm*(hnorm+lnorm)+tiny
        @check e<1//2
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
    @check length(E)==sum(dims)==1<<L
    (;E=sort!(E),dimensions=dims,blocks=records,seconds=time()-start)
end
"Floating log partition function of a spectrum (numerical scan only, not a proof input)."
logz(E,b)=let c=minimum(E); -b*c+log(sum(exp.(-b.*(E.-c)))) end
"Floating defects (1-Pₛ(n,b), 1-P_τ(2n,b)) from numerical spectra; used to rank seeds."
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
"""
    seed_enclosure(small, large, b)

Exact rational upper bounds `es ≥ 1-Pₛ(n,b)` and `et ≥ 1-P_τ(2n,b)` from the enclosed
spectra of the rings of n and 2n sites, with the three log-partition intervals.
"""
function seed_enclosure(small,large,b)
    a=enclosed_logz(small,b); c=enclosed_logz(large,b); d=enclosed_logz(large,2b)
    es=exact(sup(1-exp(c-2a))); et=exact(sup(1-exp(d-2c)))
    @check 0<=es<1 && 0<=et<1
    (;es,et,logs=[a,c,d])
end
"Paths of the 30 per-point certificates (gzip-compressed JSON, `results/g*.json.gz`)."
result_paths()=sort(filter(p->startswith(basename(p),"g") && endswith(p,".json.gz"),
    readdir(joinpath(ROOT,"results");join=true)))
"Read a per-point certificate; `.gz` files are decompressed with the system `gzip`."
read_point(path)=JSON3.read(endswith(path,".gz") ? read(`gzip -dc $path`,String) :
    read(path,String),Dict{String,Any})
"""
    hashes()

SHA-256 of this model's Julia sources (top level, `src`, `scripts`, `test`) and of the core
`julia/src`, `julia/Project.toml` and `julia/Manifest.toml`.
"""
function hashes()
    own=Tuple(filter(p->endswith(p,".jl"),readdir(ROOT)))
    merge(source_hashes(ROOT;dirs=("src","scripts","test"),files=own),
        Dict("../../julia/"*k=>v for (k,v) in source_hashes(joinpath(ROOT,"../../julia");
        dirs=("src",),files=("Project.toml","Manifest.toml"))))
end
"Run metadata (versions, threads, source hashes) for a result written by these scripts."
function metadata(config;certified=false)
    run_metadata(config;source_sha256=hashes(),packages=(CertifiedGaps,IntervalArithmetic,JSON3),
        status=certified ? "certified" : "numerical, not certified",certified)
end

"Exact seed rows (g, h, n₀, β₀, εₛ, εₜ, bound) read from the per-point certificates."
function seed_rows()
    rows=Dict{String,Any}[]
    for path in result_paths()
        d=read_point(path); c=d["certificate"]
        push!(rows,Dict{String,Any}("g"=>d["g"],"h"=>d["h"],"n0"=>c["n"],"beta0"=>c["beta"],
            "epsilon_s"=>c["epsilon_s"],"epsilon_t"=>c["epsilon_t"],
            "bound_exact_binary"=>c["bound_exact_binary"],"file"=>basename(path)))
    end
    sort!(rows;by=r->(rat(r["g"]),rat(r["h"])))
end

"""
    write_seed_table()

Write `results/seeds.json`, the exact seed data of all certificates in one small file, so that
downstream users (e.g. the Lean statement) need not decompress the spectra.
"""
function write_seed_table()
    write_result(joinpath(ROOT,"results","seeds.json"),Dict(
        "hamiltonian"=>"H_L = -sum_j Z_j Z_{j+1} - g sum_j X_j - h sum_j Z_j (Pauli, periodic)",
        "rp"=>"-h_bond = Z (x) Z with C = Z real; onsite o = -g X - h Z real symmetric",
        "epsilon_s"=>"upper bound on 1 - Z_{2n0}(beta0)/Z_{n0}(beta0)^2",
        "epsilon_t"=>"upper bound on 1 - Z_{2n0}(2 beta0)/Z_{2n0}(beta0)^2",
        "validity"=>"every even L >= n0 sites","points"=>seed_rows()))
end
end
