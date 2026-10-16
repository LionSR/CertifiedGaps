# Exact audit of internal reflection, alternating orientations, and larger blocks.
module ReflectionConventions
using LinearAlgebra
using CertifiedGaps.Spins: exact_spin_half
export audit_conventions
const Q=Rational{BigInt}
function audit_conventions()
    I2=Matrix{Q}(I,2,2)
    # Real event channels for the sublattice-rotated elementary bond.
    S=exact_spin_half(Q); channels=[S.x,S.iy,S.z]
    bond=-sum(kron(C,C) for C in channels)
    results=Any[]
    for c in (2,4)
        d=2^c; ident=Matrix{Q}(I,d,d)
        localop(C,j)=reduce(kron,[k==j ? C : I2 for k in 1:c])
        A=[localop(C,1) for C in channels]
        B=[localop(C,c) for C in channels]
        h=-sum(kron(B[k],A[k]) for k in 1:3)
        swap=[d*(b-1)+a for a in 1:d for b in 1:d]
        defect=sum(abs2,h-h[swap,swap])
        @assert defect==3*d^2//8 # Unit boundary coupling; multiply by w^2.
        reversebits(s)=sum(((s>>(j-1))&1)<<(c-j) for j in 1:c)
        reflect=[reversebits(s)+1 for s in 0:d-1]
        @assert reflect[reflect]==collect(1:d)
        @assert all(B[k][reflect,reflect]==A[k] for k in 1:3)
        for delta in (big(1)//2,big(3)//10,big(1)//4,big(1)//5,big(3)//20,big(1)//10)
            D=zeros(Q,d,d)
            for j in 1:c-1
                weight=1+(isodd(j) ? delta : -delta)
                D-=weight*sum(localop(C,j)*localop(C,j+1) for C in channels)
            end
            @assert D==D[reflect,reflect]
            @assert D==transpose(D)
        end
        if c==2
            F=ident[:,reflect]
            right=kron(ident,F); left=kron(F,ident)
            @assert right'*h*right == -sum(kron(C,C) for C in B)
            @assert left'*h*left == -sum(kron(C,C) for C in A)
            @assert left'*h*left != right'*h*right
        end
        push!(results,(;block_spins=c,local_dimension=d,unit_boundary_swap_defect=defect,
            reflected_form=true,uniform_identical_form=false))
    end
    results
end
end
