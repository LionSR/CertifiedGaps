# Untrusted exact-integer generator and independent unpacked cross-check.
# No Julia result is a Lean theorem. Float renderings are NON-CERTIFYING.
using LinearAlgebra
const ROOT = dirname(@__DIR__)
repunit(w,d)=((big(1)<<(w*d))-1)÷((big(1)<<w)-1)
pack(w,x)=sum(BigInt(x[i])<<(w*(i-1)) for i in eachindex(x);init=big(0))
lane(w,i,x)=(x>>(w*i))&((big(1)<<w)-1)
function mask(w,n,j,t)
    p=3^j
    (((big(1)<<(w*3^n))-1)÷((big(1)<<(w*3*p))-1)*((big(1)<<(w*p))-1))<<(w*t*p)
end
move(w,a,b,x)=a>=b ? x<<(w*(a-b)) : x>>(w*(b-a))
function applypacked(w,n,x)
    y=big(0)
    for j in 0:n-1
        k=mod(j+1,n)
        j0,j1,j2=[mask(w,n,j,t) for t in 0:2]
        k0,k1,k2=[mask(w,n,k,t) for t in 0:2]
        y+=(x&(j1|k1))+2*(x&((j0&k2)|(j2&k0)))
        y+=move(w,3^j,3^k,x&((j0|j1)&(k1|k2)))
        y+=move(w,3^k,3^j,x&((j1|j2)&(k0|k1)))
    end
    y
end
function reversepacked(w,d,x)
    levels=ndigits(d-1,base=2); size=1<<levels
    for level in 0:levels-1
        g=1<<level
        lo=(((big(1)<<(w*size))-1)÷((big(1)<<(w*2g))-1))*((big(1)<<(w*g))-1)
        x=((x&lo)<<(w*g))|((x>>(w*g))&lo)
    end
    x>>(w*(size-d))
end
dotpacked(w,d,x,y)=lane(w,d-1,x*reversepacked(w,d,y))
neel(n)=sum(2*3^j for j in 0:n-1 if isodd(j))

function sector(n)
    basis=[s for s in 0:3^n-1 if sum(digits(s,base=3,pad=n))==n]
    lookup=Dict(s=>i for (i,s) in enumerate(basis))
    rows=Vector{Vector{Pair{Int,Int}}}()
    for s in basis
        ds=digits(s,base=3,pad=n)
        row=Pair{Int,Int}[]; diag=n
        for j in 0:n-1
            k=mod(j+1,n)
            diag-=(ds[j+1]-1)*(ds[k+1]-1)
            for δ in (-1,1)
                0<=ds[j+1]+δ<=2 && 0<=ds[k+1]-δ<=2 || continue
                push!(row,lookup[s+δ*(3^j-3^k)]=>1)
            end
        end
        push!(row,lookup[s]=>diag); push!(rows,row)
    end
    basis,lookup,rows
end
applyrows(rows,v)=[sum(a*v[j] for (j,a) in row) for row in rows]

"""
    generate(n, q; write=true)

Power-filter trial (nI - H_M)^q e_Neel on the even untwisted n-ring: compute the exact
moments G, A, C twice (sparse rows and packed lanes), derive the Ritz/residual witnesses, and
(if `write`) emit `src/SeedPrototype/H<n>.lean`. Returns the witness row and the Lean text.
"""
function generate(n,q; write=true)
    basis,lookup,rows=sector(n)
    v=zeros(BigInt,length(basis)); v[lookup[neel(n)]]=1
    for _ in 1:q
        v=applyrows(rows,v)
    end
    av=applyrows(rows,v)
    g=dot(v,v); a=dot(v,av); c=dot(av,av)
    width=2*(q+1)*ndigits(4n-1,base=2)+ndigits(3^n-1,base=2)+2
    # Uniform row-sum bound controls ALL convolution coefficients, not just the center.
    @assert big(3)^n*big(4n)^(2*(q+1)) < big(2)^width
    pv=big(1)<<(width*neel(n))
    for _ in 1:q
        pv=applypacked(width,n,pv)
    end
    pav=applypacked(width,n,pv)
    @assert all(lane(width,s,pv)==v[i] for (i,s) in enumerate(basis))
    @assert dotpacked(width,3^n,pv,pv)==g
    @assert dotpacked(width,3^n,pv,pav)==a
    @assert dotpacked(width,3^n,pav,pav)==c
    den=10^6
    num=fld(den*(a-n*g),g)
    epsden=10^6
    epsnum=isqrt(cld(epsden^2*(c*g-a*a),g*g))+1
    name="H$(n)"
    text = """
        import SeedPrototype.Packed
        set_option maxRecDepth 100000
        set_option maxHeartbeats 0
        namespace SeedPrototype
        -- Arithmetic witnesses, NOT a theorem about ordered physical eigenvalues.
        theorem h$(n)_ritz_witness : ritzCheck $width $n $q $num $den $epsnum $epsden = true := by
          decide +kernel
        theorem h$(n)_lane_bound : $((3^n)) * (($(4n) : Nat)^$(2*(q+1))) < 2^$width := by
          decide +kernel
        #print axioms h$(n)_ritz_witness
        end SeedPrototype
        """
    write && Base.write(joinpath(ROOT,"src","SeedPrototype","$name.lean"),text)
    return (n,q,width,num,den,epsnum,epsden,g,a,c),text
end

"""(n, filter degree) pairs of the kernel-checked Ritz witnesses."""
const WITNESS_SIZES = ((6,24),(8,32),(10,16))

if abspath(PROGRAM_FILE) == @__FILE__
open(joinpath(ROOT,"results","exact_witnesses.tsv"),"w") do io
    println(io,"n\tdegree\tlane_bits\tritz_negative_num\tritz_den\tresidual_num\tresidual_den",
        "\tG\tA_filter\tC_filter")
    for (n,q) in WITNESS_SIZES
        row,_=generate(n,q)
        println(io,join(row,'\t')); flush(io)
        println("Generated H",n," degree ",q," lanes ",row[3]," bits")
    end
end
end
