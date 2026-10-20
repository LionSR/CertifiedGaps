# Float64 results in this script are NON-CERTIFYING.
using LinearAlgebra, SparseArrays, Printf
BLAS.set_num_threads(1)
const OUT = joinpath(@__DIR__, "..", "..", "results", "spin1_kernel")

function spinmatrix(n; mag=nothing, periodic=true, twist=false)
    basis = [s for s in 0:3^n-1 if mag === nothing || sum(digits(s, base=3, pad=n)) == n+mag]
    lookup = Dict(s=>i for (i,s) in enumerate(basis))
    ii, jj, vv = Int[], Int[], Float64[]
    shifts = Set{Int}()
    bondshifts = [Set{Int}() for _ in 1:n]
    for (col,s) in enumerate(basis)
        ds = digits(s, base=3, pad=n)
        diagonal = 0
        for j in 0:(periodic ? n-1 : n-2)
            k = mod(j+1,n)
            diagonal += (ds[j+1]-1)*(ds[k+1]-1)
            for delta in (-1,1)
                0 <= ds[j+1]+delta <= 2 || continue
                0 <= ds[k+1]-delta <= 2 || continue
                t = s+delta*(3^j-3^k)
                row = lookup[t]
                push!(ii,row); push!(jj,col)
                push!(vv, twist && j==n-1 ? -1 : 1)
                push!(shifts,row-col); push!(bondshifts[j+1],row-col)
            end
        end
        push!(ii,col); push!(jj,col); push!(vv,diagonal)
    end
    return sparse(ii,jj,vv,length(basis),length(basis)), basis, shifts, bondshifts
end

open(joinpath(OUT,"structure.tsv"),"w") do io
    println(io,"n\tfull_dimension\tsz0_dimension\tcompact_distinct_shifts\tcompact_per_bond",
        "\tfull_lane_MiB_at_128_bits")
    for n in (6,8,10,12)
        H, basis, shifts, bs = spinmatrix(n; mag=0)
        @printf(io,"%d\t%d\t%d\t%d\t%s\t%.4f\n",n,3^n,length(basis),length(shifts),
            join(length.(bs),","),16*3^n/2^20)
    end
end

open(joinpath(OUT,"anchor_noncertifying.tsv"),"w") do io
    println(io,"# NON-CERTIFYING Float64; H0 = two open half chains minus 4I")
    println(io,"n\topen_half_E0\tH0_E0\tshifted_cutoff\tcount_H0_below_cutoff",
        "\tcount_H0_below_minus_na")
    a=700741/500000; beta=49/4; tol=3e-5
    for n in (6,8,10,12)
        A,_,_,_ = spinmatrix(n÷2; periodic=false)
        ev=eigvals(Symmetric(Matrix(A)))
        rho= -n*a + log(3^n/tol)/beta
        pair=[x+y-4 for x in ev for y in ev]
        @printf(io,"%d\t%.12f\t%.12f\t%.12f\t%d\t%d\n",n,ev[1],2ev[1]-4,rho,
            count(<(rho),pair),count(<(-n*a),pair))
    end
end

open(joinpath(OUT,"levels_noncertifying.tsv"),"w") do io
    println(io,"# NON-CERTIFYING Float64; lowest six Sz=0 eigenvalues")
    for n in (6,8), twist in (false,true)
        A,_,_,_=spinmatrix(n;mag=0,twist)
        ev=eigvals(Symmetric(Matrix(A)))
        println(io,n,'\t',twist,'\t',join(ev[1:6],"\t"))
    end
end

open(joinpath(OUT,"anchor_separation_noncertifying.tsv"),"w") do io
    println(io,"# NON-CERTIFYING Float64; ideal Ritz values equal exact eigenvalues")
    println(io,"n\td\tmax_over_k_H0_lambda_kplus1_minus_H_lambda_k\tbest_k\tviable_k_count")
    for n in (6,8)
        ev=Float64[]
        for mag in -n:n
            H,_,_,_=spinmatrix(n;mag)
            append!(ev,eigvals(Symmetric(Matrix(H))))
        end
        sort!(ev)
        A,_,_,_=spinmatrix(n÷2;periodic=false)
        e0=eigvals(Symmetric(Matrix(A)))
        anchor=sort([x+y-4 for x in e0 for y in e0])
        margins=anchor[2:end]-ev[1:end-1]
        @printf(io,"%d\t%d\t%.12f\t%d\t%d\n",n,length(ev),maximum(margins),argmax(margins),
            count(>(1e-9),margins))
    end
end
