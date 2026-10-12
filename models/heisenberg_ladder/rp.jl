using LinearAlgebra

"Verify the checkerboard rotation and cell RP decomposition over Rational{BigInt}."
function exact_rp(J::Rational)
    T = Rational{BigInt}
    x = T[0 1; 1 0] / 2
    y = T[0 1; -1 0] / 2 # i Sy, real antisymmetric
    z = T[1 0; 0 -1] / 2
    id = Matrix{T}(I, 2, 2)
    u = 2y # exp(i pi Sy), exactly orthogonal
    h = kron(x,x) - kron(y,y) + kron(z,z)
    hr = -kron(x,x) - kron(y,y) - kron(z,z)
    @assert kron(u,id)*h*transpose(kron(u,id)) == hr
    @assert kron(id,u)*h*transpose(kron(id,u)) == hr
    C = [kron(a,id) for a in (x,y,z)]
    append!(C, [kron(id,a) for a in (x,y,z)])
    # Bond between two one-rung cells; each of the two legs has coefficient 1.
    hb = -sum(kron(c,c) for c in C)
    original = sum(kron(kron(a,id),kron(a,id)) +
                   kron(kron(id,a),kron(id,a)) for a in (x,z)) -
               kron(kron(y,id),kron(y,id)) - kron(kron(id,y),kron(id,y))
    ua, ub = kron(u,id), kron(id,u)
    @assert kron(ua,ub)*original*transpose(kron(ua,ub)) == hb
    @assert kron(ub,ua)*original*transpose(kron(ub,ua)) == hb
    onsite = J*hr
    @assert onsite == transpose(onsite)
    @assert ua*(J*h)*transpose(ua) == onsite == ub*(J*h)*transpose(ub)
    @assert -hb == sum(kron(c,c) for c in C)
    return true
end
