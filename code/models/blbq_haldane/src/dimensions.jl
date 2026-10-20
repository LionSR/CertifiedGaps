"Exact dimension of the spin-1 magnetization-M space."
function charge_count(L,M)
    abs(M)>L && return 0
    coeff=[1]
    for _ in 1:L
        next=zeros(Int,length(coeff)+2)
        for i in eachindex(coeff), j in 0:2; next[i+j]+=coeff[i]; end
        coeff=next
    end
    coeff[M+L+1]
end

"Möbius function μ(n)."
function mobius(n)
    sign=1
    for p in 2:n
        n%p==0 || continue
        n÷=p;sign=-sign
        n%p==0 && return 0
    end
    sign
end

"Exact translation-sector dimension from integer Ramanujan sums."
function dimension(L,M,k)
    value=0
    for r in 1:L
        L%r==0 && M%r==0 || continue
        ramanujan=sum(d*mobius(r÷d) for d in 1:r if gcd(r,k)%d==0)
        value+=ramanujan*charge_count(L÷r,M÷r)
    end
    @assert value%L==0
    value÷L
end
