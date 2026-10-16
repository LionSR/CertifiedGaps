# Exact reflection positivity of the checkerboard-rotated ladder. Spectra are computed in the
# original SU(2) basis; these identities are what make that basis legitimate.

"""
    exact_rp(J) -> true

Verify over `Rational{BigInt}` that the checkerboard rotation u = exp(iπSʸ) = 2iSʸ maps the
original rung and leg bonds to the reflection-positive cell model, for both alternating
orientations, with `-h_cell = Σ_C C ⊗ C` for the six real factors
`C ∈ {X⊗I, A⊗I, Z⊗I, I⊗X, I⊗A, I⊗Z}` (X = Sˣ, A = iSʸ, Z = Sᶻ) and a real symmetric
on-rung term. Throws on any failed identity.
"""
function exact_rp(J::Rational)
    S = exact_spin_half(Q)
    x, y, z, id = S.x, S.iy, S.z, S.id # y = i Sʸ, real antisymmetric
    u = 2y # exp(iπSʸ), exactly orthogonal
    # S·S written with real matrices: Sʸ⊗Sʸ = -(iSʸ)⊗(iSʸ).
    h = kron(x, x) - kron(y, y) + kron(z, z)
    hr = -kron(x, x) - kron(y, y) - kron(z, z)
    @assert kron(u, id) * h * transpose(kron(u, id)) == hr
    @assert kron(id, u) * h * transpose(kron(id, u)) == hr
    C = [kron(a, id) for a in (x, y, z)]
    append!(C, [kron(id, a) for a in (x, y, z)])
    # Bond between two one-rung cells; each of the two legs has coefficient 1.
    hb = -sum(kron(c, c) for c in C)
    original = sum(kron(kron(a, id), kron(a, id)) +
                   kron(kron(id, a), kron(id, a)) for a in (x, z)) -
               kron(kron(y, id), kron(y, id)) - kron(kron(id, y), kron(id, y))
    ua, ub = kron(u, id), kron(id, u)
    # Alternating rungs are rotated on opposite legs, so both bond orientations occur.
    @assert kron(ua, ub) * original * transpose(kron(ua, ub)) == hb
    @assert kron(ub, ua) * original * transpose(kron(ub, ua)) == hb
    onsite = J * hr
    @assert onsite == transpose(onsite)
    @assert ua * (J * h) * transpose(ua) == onsite == ub * (J * h) * transpose(ub)
    @assert -hb == sum(kron(c, c) for c in C)
    return true
end

"""
    rp_model(J) -> (model, rp)

The rotated one-rung cell model as a core `RingModel` (d = 4) together with its explicit
real reflection-positive factorisation, checked by `CertifiedGaps.RingModels.verify_rp`.
"""
function rp_model(J)
    S = exact_spin_half(Q)
    x, a, z, id = S.x, S.iy, S.z, S.id
    factors = vcat([kron(c, id) for c in (x, a, z)], [kron(id, c) for c in (x, a, z)])
    bond = -sum(kron(c, c) for c in factors)
    onsite = -J * (kron(x, x) + kron(a, a) + kron(z, z))
    model = RingModel(4, bond, onsite; name = "checkerboard-rotated Heisenberg ladder")
    return model, rp_certificate(model; factors)
end
