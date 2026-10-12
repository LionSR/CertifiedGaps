# Independent exact-arithmetic check of floating residual allowances on all L=6 blocks.
include("certify.jl")
using .LongitudinalCertificate, LinearAlgebra, CertifiedGaps.Results
const C = LongitudinalCertificate
q(x) = Rational{BigInt}(x)
cq(z) = complex(q(real(z)), q(imag(z)))
start = time()
count = 0
largest_e_ratio = big(0)//big(1)
largest_r_ratio = big(0)//big(1)
for g in C.GS, h in C.HS, k in 0:5
    H, a = C.block(C.orbits(6), k, g, h; rigorous=true)
    F = eigen(Hermitian(H))
    Q, A, lam = cq.(F.vectors), cq.(H), q.(F.values)
    # These matrix products, squared norms, and comparisons are entirely exact.
    R, E = A*Q-Q.*transpose(lam), Q'*Q-I
    r2, e2 = sum(abs2, R), sum(abs2, E)
    enc = C.enclosed_eigen(H, a)
    @assert r2 <= q(enc.residual)^2
    @assert e2 <= q(enc.orthogonality)^2
    global largest_e_ratio = max(largest_e_ratio, e2/q(enc.orthogonality)^2)
    global largest_r_ratio = max(largest_r_ratio, r2/q(enc.residual)^2)
    global count += 1
end
@assert count == 180
elapsed = time()-start
write_result(joinpath(C.ROOT, "results", "independent_residual_audit.json"), Dict(
    "status" => "passed", "ring_length" => 6, "parameter_points" => 30,
    "blocks" => count, "exact_rational_inequalities" => 2count,
    "largest_squared_residual_ratio" => C.str(largest_r_ratio),
    "largest_squared_orthogonality_ratio" => C.str(largest_e_ratio),
    "seconds" => elapsed,
    "scope" => "Exact rational checks of implemented error allowances; not a complete spectral replay",
    "metadata" => C.metadata(Dict("kind" => "independent exact rational residual audit"))
))
println("Passed ", 2count, " exact rational inequalities in ", elapsed, " seconds")
