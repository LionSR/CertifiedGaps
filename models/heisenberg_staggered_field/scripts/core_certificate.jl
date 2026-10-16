# Independent certificate of the h = 1 seed with the shared core (CertifiedGaps.GapCertificates):
# its own exact RP identity, verified block spectra and exact rational chain, in the rotated
# frame where only S^z parity is conserved. Produces results/core_certificate_1_1.json.
#     julia --project=models/heisenberg_staggered_field \
#         models/heisenberg_staggered_field/scripts/core_certificate.jl        (≈ 12 min)
# test/runtests.jl checks its exact parts; `check.jl --full` also re-verifies its spectra.

using CertifiedGaps.RingModels, CertifiedGaps.GapCertificates, CertifiedGaps.Results
using CertifiedGaps.Spins: exact_spin_half

# Spin-½ operators with the real antisymmetric J = iSʸ.
const X, J, Z = let S = exact_spin_half(Rational{Int}); (S.x, S.iy, S.z) end

"Rotated-frame staggered Heisenberg ring: h = −(X⊗X + J⊗J + Z⊗Z), o = −h_s Z."
staggered(hs) = RingModel(2, -(kron(X, X) + kron(J, J) + kron(Z, Z)), -hs * Z;
    charge = [0, 1], modulus = 2, name = "spin-1/2 Heisenberg, staggered field $hs, rotated")

if abspath(PROGRAM_FILE) == @__FILE__
    model = staggered(1)
    cert = certify_gap(model; n0 = 8, β0 = 21 // 4,
        rp = rp_certificate(model; factors = [X, J, Z]), progress = true)
    path = write_result(joinpath(dirname(@__DIR__), "results", "core_certificate_1_1.json"), cert)
    println(path, ": gap > ", cert["gap_lower_bound_decimal"])
end
