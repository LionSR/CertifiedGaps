"""
    CertifiedGaps.CertifiedMPO.AdjointTrace

Conditional trace enclosure for a sequence of approximate Hermitian squares or reflected
Gram steps (`trace_interval`). Moved verbatim from spin2 `Spin2AdjointTrace`.
"""
module AdjointTrace

using ..CertifiedMPO: F

"""
Conditional trace enclosure for a sequence of approximate Hermitian squares.

Let A[j+1]=A[j]^2/kappa[j], with Hermitian A and B, and G[end]=B[end].
Define R[j]=B[j]^2/kappa[j]-B[j+1] and
T[j]=(B[j]*G[j+1]+G[j+1]*B[j])/kappa[j]-G[j], with Hermitian G.
The signed center is
  norm(B[end],HS)^2 + 2inner(G[1],A[1]-B[1]) + 2sum(inner(G[j+1],R[j])).
Inputs enclose this center and bound the HS primal errors, HS dual defects,
and operator norms of G[j+1]. All these inputs require separate certification.

Telescoping gives an uncertain remainder bounded by
  2sum(tau[j]*epsilon[j]+gamma[j]*epsilon[j]^2/kappa[j]),
plus the nonnegative terminal squared error, at most epsilon[end]^2.

The same scalar enclosure applies to non-Hermitian candidates and reflected
Gram steps F_j(A)=S_j(A'*A)/kappa[j], where S_j is identity or a self-adjoint
unitary spatial reflection. Use the real HS inner product and replace R[j]
by F_j(B[j])-B[j+1] and T[j] by
  B[j]*(S_j(G[j+1])+S_j(G[j+1])')/kappa[j]-G[j].
The nonlinear remainder is S_j(D[j]'*D[j])/kappa[j]; its pairing with
G[j+1] has absolute value at most norm(G[j+1],op)*norm(D[j],HS)^2/kappa[j].
"""
function trace_interval(
    center_lower::F,
    center_upper::F;
    errors,
    dual_defects,
    dual_operator_bounds,
    normalizations,
    initial_uncertainty = F(0),
)
    steps = length(normalizations)
    @assert length(errors) == steps+1
    @assert length(dual_defects) == length(dual_operator_bounds) == steps
    @assert center_lower <= center_upper && initial_uncertainty >= 0
    @assert all(>=(0), errors) && all(>=(0), dual_defects) && all(>=(0), dual_operator_bounds)
    @assert all(>(0), normalizations)
    cross = 2sum(dual_defects[j]*errors[j] for j in 1:steps; init = F(0))
    quadratic =
        2sum(dual_operator_bounds[j]*errors[j]^2/normalizations[j] for j in 1:steps; init = F(0))
    radius = initial_uncertainty+cross+quadratic
    return (;
        lower = center_lower-radius,
        upper = center_upper+radius+last(errors)^2,
        cross,
        quadratic,
        terminal = last(errors)^2,
        radius,
    )
end

end
