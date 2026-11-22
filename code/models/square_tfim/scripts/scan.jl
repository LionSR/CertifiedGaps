# Floating-point diagnostics only. Certificates are produced by certify.jl.
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
using CertifiedGaps.Results: write_result
using LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "SquareIsing.jl"))
using .SquareIsing

BLAS.set_num_threads(1)
g = isempty(ARGS) ? 24 : parse(Int, ARGS[1])
f(z) = z^2 / (2(1-z)^2)
function step(e)
    x, y, t = e
    a = f(1 - (1-x)^4 * (1-y)^2 * (1-t))
    a < 1 || return (Inf, Inf, Inf)
    b = f(1 - (1-y)^4 * (1-t)^2 * (1-a))
    b < 1 || return (Inf, Inf, Inf)
    c = f(1 - (1-t)^4 * (1-a)^2 * (1-b))
    return (a, b, c)
end
spectra = []
for (n, m) in ((2, 2), (4, 2), (4, 4))
    elapsed = @elapsed s = float_ring_spectrum(row_model(m, g), n)
    push!(spectra, s.E)
    println("g=$g torus=$n,$m dimension=$(length(s.E)) seconds=$elapsed")
    flush(stdout)
end
logz(e, b) = -b * first(e) + log(sum(exp(-b * (x-first(e))) for x in e))
rows = []
for b in 2.0:0.125:12.0
    beta = b/g
    a, c, d = [logz(e, beta) for e in spectra]
    e = (-expm1(c-2a), -expm1(d-2c), -expm1(logz(spectra[3], 2beta)-2d))
    q = (1-e[1])^4 * (1-e[2])^2 * (1-e[3])
    bound = q > 0.5 ? log(q/(1-q))/beta : 0.0
    current = e
    steps = 0
    while maximum(current) > 0.001 && maximum(current) < 1 && steps < 12
        current = step(current)
        steps += 1
        q = prod((1-current[i])^([4,2,1][i]) for i in 1:3)
        q > 0.5 && (bound = min(bound, log(q/(1-q))/(2^steps*beta)))
    end
    success = maximum(current) <= 0.001
    success && (bound = min(bound, log(1/(25maximum(current)))/(2^steps*beta)))
    push!(rows, (; b, beta, defects = e, success, steps, bound = success ? bound : 0.0))
end
best = rows[argmax([r.bound for r in rows])]
println("BEST: $best")
write_result(joinpath(@__DIR__, "..", "results", "scan_g$g.json"),
    Dict("status" => "floating-point diagnostics, not certified", "g" => string(g),
        "rows" => rows, "best" => best))
