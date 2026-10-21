"""
    Spin1Kernel

Study `spin1_kernel`: packed-`Nat` arithmetic for the spin-1 seed in the Lean kernel. Untrusted
exact-integer generator of the Ritz/residual witnesses checked by `lean/SeedKernel/Spin1/H*.lean`
(with an independent unpacked cross-check), the witness inventory `results/spin1_kernel/
witnesses.json`, and its checker. No Julia result is a Lean theorem; float renderings are
NON-CERTIFYING.
"""
module Spin1Kernel

using LinearAlgebra
using CertifiedGaps.Results: write_result, sha256_file, run_metadata, read_exact_json
using ..Common: MODEL, LEAN, results_dir, study_hashes, export_all

const OUT = results_dir("spin1_kernel")
"Directory of the generated `H<n>.lean` witness files."
const LEAN_SPIN1 = joinpath(LEAN, "SeedKernel", "Spin1")

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
    lean_text(n, q, width, num, den, epsnum, epsden)

The Lean source of `SeedKernel/Spin1/H<n>.lean`: the kernel-checked `ritzCheck` witness and the
lane-width bound, as a function of the witness row alone.
"""
function lean_text(n,q,width,num,den,epsnum,epsden)
    """
    import SeedKernel.Spin1.Packed
    set_option maxRecDepth 100000
    set_option maxHeartbeats 0
    namespace SeedKernel.Spin1
    -- Arithmetic witnesses, NOT a theorem about ordered physical eigenvalues.
    theorem h$(n)_ritz_witness : ritzCheck $width $n $q $num $den $epsnum $epsden = true := by
      decide +kernel
    theorem h$(n)_lane_bound : $((3^n)) * (($(4n) : Nat)^$(2*(q+1))) < 2^$width := by
      decide +kernel
    #print axioms h$(n)_ritz_witness
    end SeedKernel.Spin1
    """
end

"""
    generate(n, q; write=true)

Power-filter trial (nI - H_M)^q e_Neel on the even untwisted n-ring: compute the exact
moments G, A, C twice (sparse rows and packed lanes), derive the Ritz/residual witnesses, and
(if `write`) emit `lean/SeedKernel/Spin1/H<n>.lean`. Returns the witness row and the Lean text.
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
    text=lean_text(n,q,width,num,den,epsnum,epsden)
    write && Base.write(joinpath(LEAN_SPIN1,"H$n.lean"),text)
    return (n,q,width,num,den,epsnum,epsden,g,a,c),text
end

"""(n, filter degree) pairs of the kernel-checked Ritz witnesses."""
const WITNESS_SIZES = ((6,24),(8,32),(10,16))

"""
    lean_sources()

The Lean and Lake files the spin-1 witnesses depend on (shared lanes, `Spin1/`, benchmarks,
audit, Lake configuration), relative to MODEL.
"""
function lean_sources()
    files = String[]
    seed = joinpath(LEAN, "SeedKernel")
    for name in ("Packed", "PackingProofs", "Obligations", "PackedTests", "Audit")
        push!(files, joinpath(seed, name * ".lean"))
    end
    for dir in ("Spin1", "Bench"), name in readdir(joinpath(seed, dir))
        endswith(name, ".lean") && push!(files, joinpath(seed, dir, name))
    end
    append!(files, joinpath.(LEAN, ["SeedKernel.lean", "lakefile.toml", "lean-toolchain",
        "lake-manifest.json"]))
    return sort([relpath(f, MODEL) for f in files])
end

"Hashes of everything the checker trusts or replays: Lean sources and this study's Julia."
function hashed_sources()
    hashes = study_hashes("spin1_kernel")
    for f in lean_sources()
        hashes[f] = sha256_file(joinpath(MODEL, f))
    end
    return hashes
end

"Witness rows as exact strings (rationals as `num//den`)."
witness_dict(n,q,width,num,den,epsnum,epsden,g,a,c) = Dict("n" => "$n", "degree" => "$q",
    "lane_bits" => "$width", "ritz_upper" => "-$num//$den",
    "residual_upper" => "$epsnum//$epsden", "G" => "$g", "A_filter" => "$a", "C_filter" => "$c")

"""
    certify_spin1_kernel()

Regenerate `lean/SeedKernel/Spin1/H{6,8,10}.lean`, `results/spin1_kernel/exact_witnesses.tsv`
and `results/spin1_kernel/witnesses.json` (exact witnesses, source hashes, versions). About
two minutes. This is an inventory, NOT a gap certificate.
"""
function certify_spin1_kernel()
    rows = Any[]
    open(joinpath(OUT,"exact_witnesses.tsv"),"w") do io
        println(io,"n\tdegree\tlane_bits\tritz_negative_num\tritz_den\tresidual_num",
            "\tresidual_den\tG\tA_filter\tC_filter")
        for (n,q) in WITNESS_SIZES
            row,_=generate(n,q)
            println(io,join(row,'\t')); flush(io)
            push!(rows,witness_dict(row...))
            println("Generated H",n," degree ",q," lanes ",row[3]," bits")
        end
    end
    toolchain = strip(read(joinpath(LEAN, "lean-toolchain"), String))
    status = "kernel arithmetic witnesses only; no physical eigenvalue or gap theorem"
    write_result(joinpath(OUT, "witnesses.json"), Dict(
        "status" => status,
        "witnesses" => rows,
        "run_metadata" => run_metadata(Dict("lean_toolchain" => toolchain,
            "lean_dependencies" =>
                "Std only (SeedKernel.Packed and SeedKernel.Spin1 import no Mathlib)");
            source_sha256=hashed_sources(), status, certified=false)))
end

"""
    check_spin1_kernel(; full=false)

Compare the recorded source hashes with the current files, and each `H<n>.lean` with the text
the generator emits for the stored witness row; with `full`, also recompute every witness from
scratch (about two minutes). The Lean kernel itself is run by `Common.build_seed_kernel`.
Throws on any mismatch.
"""
function check_spin1_kernel(; full=false)
    stored = read_exact_json(joinpath(OUT, "witnesses.json"))
    recorded = Dict(string(k) => string(v) for (k, v) in stored["run_metadata"]["source_sha256"])
    current = hashed_sources()
    recorded == current || error("source hashes differ from results/spin1_kernel/" *
        "witnesses.json:\n" * join(sort(collect(symdiff(Set(recorded), Set(current)))), "\n"))
    println("spin1_kernel: source hashes match ($(length(current)) files)")
    for (row, (n, q)) in zip(stored["witnesses"], WITNESS_SIZES)
        row = Dict(string(k) => string(v) for (k, v) in row)
        row["n"] == "$n" && row["degree"] == "$q" || error("unexpected witness order")
        num, den = parse.(BigInt, split(chopprefix(row["ritz_upper"], "-"), "//"))
        epsnum, epsden = parse.(BigInt, split(row["residual_upper"], "//"))
        width = parse(Int, row["lane_bits"])
        text = lean_text(n, q, width, num, den, epsnum, epsden)
        text == read(joinpath(LEAN_SPIN1, "H$n.lean"), String) ||
            error("lean/SeedKernel/Spin1/H$n.lean differs from the stored witness")
        if full
            fresh, fresh_text = generate(n, q; write=false)
            witness_dict(fresh...) == row || error("H$n witness differs from regeneration")
            fresh_text == text || error("H$n Lean text differs from regeneration")
        end
        println("spin1_kernel: H$n Ritz <= -$num/$den, residual <= $epsnum/$epsden, ",
            "$width-bit lanes", full ? " (regenerated)" : "")
    end
    return true
end

export_all(@__MODULE__)

end
