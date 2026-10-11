using LinearAlgebra
"Exact Cartesian site RP and reflection within the two-site cell."
function check_rp(delta)
    @assert -1<delta<1
    G=[BigInt[0 0 0;0 0 -1;0 1 0],BigInt[0 0 1;0 0 0;-1 0 0],BigInt[0 -1 0;1 0 0;0 0 0]]
    G=[-g for g in G]
    S=[-im*g for g in G]
    @assert sum(s*s for s in S)==2Matrix{BigInt}(I,3,3)
    @assert S[1]*S[2]-S[2]*S[1]==im*S[3]
    h=sum(kron(s,s) for s in S)
    for c in (1-delta,1+delta)
        @assert -c*h==sum(c*kron(g,g) for g in G)
    end
    P=zeros(BigInt,9,9)
    for a in 1:3,b in 1:3;P[3(b-1)+a,3(a-1)+b]=1;end
    D=(1-delta)*h
    @assert isreal(D) && D==transpose(D) && P*conj(D)*P==D
    A=[kron(Matrix{BigInt}(I,3,3),g) for g in G]
    B=[kron(g,Matrix{BigInt}(I,3,3)) for g in G]
    @assert all(P*conj(A[a])*P==B[a] for a in 1:3)
    hb=(1+delta)*sum(kron(kron(Matrix{BigInt}(I,3,3),s),kron(s,Matrix{BigInt}(I,3,3))) for s in S)
    @assert -hb==(1+delta)*sum(kron(A[a],P*conj(A[a])*P) for a in 1:3)
    true
end
