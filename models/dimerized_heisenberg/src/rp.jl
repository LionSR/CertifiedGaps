module DimerizedRP
using LinearAlgebra
export audit_rp
const Q = Rational{BigInt}
const CQ = Complex{Q}

"Exact audit of the requested real identical-factor RP form, plus the broader conjugate form."
function audit_rp(delta::Q)
    ident = Matrix{CQ}(I, 2, 2)
    sx = CQ[0 1; 1 0] / 2
    sy = CQ[0 -im; im 0] / 2
    sz = CQ[1 0; 0 -1] / 2
    spins = (sx, sy, sz)
    A = [kron(s, ident) for s in spins]
    B = [kron(ident, s) for s in spins]
    weak, strong = 1-delta, 1+delta
    bond = weak * sum(kron(B[k], A[k]) for k in 1:3)
    onsite = strong * sum(A[k]*B[k] for k in 1:3)
    swap = zeros(CQ, 16, 16)
    for a in 1:4, b in 1:4
        swap[4*(b-1)+a, 4*(a-1)+b] = 1
    end
    defect = bond - swap*bond*swap
    defect_norm2 = sum(abs2, defect)
    @assert defect_norm2 == 6*weak^2
    # A cellwise pi rotation about y on the second spin.
    J = CQ[0 1; -1 0]
    rotation = kron(ident, J)
    internal_swap = CQ[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    Drot = rotation'*onsite*rotation
    R2 = kron(rotation,rotation)
    hrot = R2'*bond*R2
    real_channels = [B[1], im*B[2], B[3]]
    theta(C) = internal_swap*conj.(C)*internal_swap
    @assert all(C -> all(iszero,imag.(C)), real_channels)
    @assert all(iszero,imag.(Drot))
    @assert internal_swap*conj.(internal_swap) == Matrix{CQ}(I,4,4)
    @assert theta(Drot) == Drot
    @assert -hrot == weak*sum(kron(C,theta(C)) for C in real_channels)
    # W = M/sqrt(2): three real triplets and an imaginary singlet.
    M = CQ[1 1 0 0; 0 0 1 im; 0 0 1 -im; 1 -1 0 0]
    @assert M'*M == 2*Matrix{CQ}(I,4,4)
    V = rotation*M # Full basis U=V/sqrt(2); evaluate rationally.
    UV = kron(V,V)
    transformed_bond = UV'*bond*UV/4
    transformed_onsite = V'*onsite*V/2
    @assert all(iszero, imag.(transformed_onsite))
    @assert transformed_onsite == transpose(transformed_onsite)
    transformed_defect = transformed_bond-swap*transformed_bond*swap
    @assert sum(abs2, transformed_defect) == defect_norm2
    # In the rotated product frame the channels are Bx, i By, Bz.
    channels = [M'*B[1]*M/2, M'*(im*B[2])*M/2, M'*B[3]*M/2]
    conjugate_rhs = weak*sum(kron(C,conj.(C)) for C in channels)
    real_rhs = weak*sum(kron(C,C) for C in channels)
    @assert -transformed_bond == conjugate_rhs
    @assert any(C -> any(!iszero, imag.(C)), channels)
    @assert -transformed_bond != real_rhs
    @assert defect_norm2 > 0
    return (; delta, weak, strong, defect_norm2, transformed_onsite, channels,
        conjugate_identity=true, reflected_real_identity=true,
        theta_invariant_onsite=true, real_identical_identity=false,
        real_rhs_error_norm2=sum(abs2,-transformed_bond-real_rhs),
        transformed_bond, basis_numerator=V, mirror_basis_numerator=M)
end
end
