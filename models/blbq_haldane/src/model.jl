using LinearAlgebra

"Exact coefficient identities, hence valid for symbolic c=cos(theta), s=sin(theta)."
function exact_rp()
    G = [Rational{BigInt}.(a) for a in (
        [0 0 0;0 0 1;0 -1 0], [0 0 -1;0 0 0;1 0 0], [0 1 0;-1 0 0;0 0 0])]
    Q = [-(G[a]*G[b]+G[b]*G[a])/2 for a in 1:3, b in 1:3]
    A = sum(kron(g,g) for g in G)
    B = sum(kron(q,q) for q in Q)
    X = -A
    @assert -X == A
    @assert -X^2 == -A/2-B
    @assert all(g -> transpose(g) == -g, G)
    @assert all(q -> transpose(q) == q, Q)
    onsite = zeros(Rational{BigInt},3,3)
    @assert onsite == transpose(onsite)
    # Independent realignment check: scalar, symmetric-traceless, skew subspaces.
    swap=zeros(Rational{BigInt},9,9)
    for a in 1:3,b in 1:3
        swap[a+3(b-1),b+3(a-1)]=1
    end
    ident=Matrix{Rational{BigInt}}(I,9,9)
    scalar=vec(Matrix{Rational{BigInt}}(I,3,3))
    @assert sum(vec(g)*vec(g)' for g in G)==ident-swap
    @assert sum(vec(q)*vec(q)' for q in Q)==(ident+swap)/2+scalar*scalar'
    return (; G,Q,X,A,B,onsite)
end

"Spin-dot-spin in the Sz basis, with exact integer entries for spin one."
function sz_bond()
    X = zeros(Int,9,9)
    for a in 0:2, b in 0:2
        j=3a+b+1
        X[j,j]=(a-1)*(b-1)
        a<2 && b>0 && (X[3(a+1)+(b-1)+1,j]=1)
        a>0 && b<2 && (X[3(a-1)+(b+1)+1,j]=1)
    end
    return X
end
