# Extra rigorous path for the twelve-spin ring. No Lanczos, truncated sectors or
# assumptions about the numerical projector eigenvectors are used.
"An exact integer multiple of a dihedral character projector on one group orbit."
function integer_projector(os,cells,k,p,reflection)
    L=2cells;flipmax=3^L-1;oi=Dict(t=>j for (j,t) in enumerate(os))
    generic=!(k==0 || 2k==cells)
    den=2*cells*2*(p==0 ? 1 : 2)
    P=zeros(Int,length(os),length(os))
    for (j,v) in enumerate(os),f in 0:(p==0 ? 0 : 1),h in 0:1
        t=f==0 ? v : flipmax-v;t=h==0 ? t : reflect_state(t,L)
        for r in 0:cells-1
            # For cells in {2,3,4,6}, twice every cosine is exactly integral.
            @assert cells in (2,3,4,6)
            angle=mod(k*r,cells)
            cosine2=angle==0 ? 2 : 2angle==cells ? -2 :
                (4angle==cells || 4angle==3cells) ? 0 :
                (6angle==cells || 6angle==5cells) ? 1 :
                (3angle==cells || 3angle==2cells) ? -1 : error("unsupported angle")
            P[oi[t],j]+=(generic ? 2 : 1)*cosine2*(f==0 ? 1 : p)*(h==0 ? 1 : reflection)
            t=translate(t,L)
        end
    end
    @assert P==P' && P*P==den*P && tr(P)%den==0
    P,den,tr(P)÷den
end
"Bound distance to an exactly orthonormal basis in the exact projector range."
function basis_distance(Pnum,den,B)
    s=size(B,1);r=size(B,2);r==0 && return 0.0
    P=Pnum/den
    R=P*B-B;E=B'*B-Matrix{Float64}(I,r,r)
    # Frobenius <= s*max-entry; all other arithmetic below is outward rounded.
    u=interval(1)/interval(big(2)^53)
    gamma=(2s+8)*u/(1-(2s+8)*u)
    # ||P_float-P_exact||_F <= s*u: exact projector entries have modulus <=1.
    b=interval(s)*maximum(abs,B)
    projection=s*interval(maximum(abs,R))+gamma*(s*b+b)+s*u*b
    orth=s*interval(maximum(abs,E))+gamma*(b*b+s)
    wdefect=orth+2b*projection+projection^2
    @assert sup(wdefect)<0.5
    dist=projection+wdefect/(1+sqrt(1-wdefect))
    sup(dist)
end
"Basis plus a rigorous distance bound to an exact orthonormal symmetry basis."
function verified_basis(states,cells,k,p,reflection)
    L=2cells;flipmax=3^L-1;idx=Dict(s=>i for (i,s) in enumerate(states))
    visited=Set{Int}();rows=Int[];cols=Int[];vals=Float64[];col=0;distance=0.0
    for s in states
        s in visited && continue
        orbit=Set{Int}()
        for f in 0:(p==0 ? 0 : 1),h in 0:1
            t=f==0 ? s : flipmax-s;t=h==0 ? t : reflect_state(t,L)
            for r in 0:cells-1;push!(orbit,t);t=translate(t,L);end
        end
        union!(visited,orbit);os=sort!(collect(orbit))
        P,den,rank=integer_projector(os,cells,k,p,reflection)
        rank==0 && continue
        F=eigen(Symmetric(P/den)); B=F.vectors[:,end-rank+1:end]
        distance=max(distance,basis_distance(P,den,B))
        for a in 1:rank
            col+=1
            for (i,t) in enumerate(os)
                push!(rows,idx[t]);push!(cols,col);push!(vals,B[i,a])
            end
        end
    end
    sparse(rows,cols,vals,length(states),col),distance
end
"Rigorous sparse projected Hamiltonian construction error, including basis error."
function projected_hamiltonian(A,den,B,distance,cells,delta)
    raw=(B'*(A*B))/den
    # Enforce exact symmetry by using the computed upper triangle in both halves.
    upper=triu(raw);H=upper+transpose(triu(raw,1))
    dim=size(B,2); maxa=maximum(diff(A.colptr)); maxb=maximum(diff(B.colptr))
    # A*B dot products have at most maxa terms, B'*AB at most maxb terms.
    # Account for both products, rational scaling and triangle copying.
    steps=2maxa+2maxb+12
    model_error=setprecision(BigFloat,128) do
        setrounding(BigFloat,RoundUp) do
            gamma=BigFloat(big(steps)//(big(2)^53-steps))
            rounderror=2gamma*frob_up(nonzeros(A))*frob_up(nonzeros(B))^2/den
            # Sum_j ||J_j S.S|| = 2 sum_j J_j = 4*cells.
            d=BigFloat(distance)
            rounderror+4cells*(2d+d*d)
        end
    end
    H,model_error
end
"Complete verified real symmetry sectors for a large ring."
function verified_block_spectrum(cells,delta)
    groups=charge_groups(cells); vals=Float64[]; radii=Rational{BigInt}[];records=Any[]
    start=time()
    for m in 0:2cells
        states=groups[m+1];A,den=integer_hamiltonian(cells,delta,states)
        for k in 0:fld(cells,2),p in (m==0 ? (-1,1) : (0,)),reflection in (k==0 || 2k==cells ? (-1,1) : (1,))
            begin_time=time()
            B,dist=verified_basis(states,cells,k,p,reflection);dim=size(B,2)
            dim==0 && continue
            H,model_error=projected_hamiltonian(A,den,B,dist,cells,delta)
            lam,radius,record=eigen_enclosure(H,1)
            error_rational=setrounding(BigFloat,RoundUp) do
                ceil(BigInt,model_error*BigFloat(10)^12)//big(10)^12
            end
            radius+=error_rational
            mult=(m==0 ? 1 : 2)*(k==0 || 2k==cells ? 1 : 2)
            record["charge"]=m;record["momentum"]=k;record["spin_flip"]=p
            record["reflection"]=reflection;record["multiplicity"]=mult
            record["basis_distance_upper"]=ratstr(Rational{BigInt}(dist))
            record["model_error_upper"]=ratstr(error_rational)
            record["eigenvalue_radius"]=ratstr(radius);record["seconds"]=time()-begin_time
            push!(records,record)
            for _ in 1:mult;append!(vals,lam);append!(radii,fill(radius,length(lam)));end
            println("VERIFIED cells=$cells M=$m k=$k p=$p R=$reflection dim=$dim radius=$radius seconds=$(record["seconds"])");flush(stdout)
        end
    end
    @assert length(vals)==9^cells==sum(r["dimension"]*r["multiplicity"] for r in records)
    vals,radii,records,time()-start
end
