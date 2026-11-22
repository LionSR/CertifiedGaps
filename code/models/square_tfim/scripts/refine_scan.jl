# Floating-point temperature search; certify.jl independently encloses the chosen points.
using LinearAlgebra
using CertifiedGaps.Results
include(joinpath(@__DIR__, "..", "src", "Certificate.jl"))
using .Certificate

BLAS.set_num_threads(1)
root = dirname(@__DIR__)

function approximate_logz(spectrum, beta)
    e0 = minimum(minimum(block.centers) for block in spectrum)
    return -beta*e0 + log(sum(block.multiplicity *
        sum(exp(-beta*(e-e0)) for e in block.centers) for block in spectrum))
end

function approximate_bound(defects, beta)
    f(z) = z^2 / (2*(1-z)^2)
    zq(e) = -expm1(4log1p(-e[1]) + 2log1p(-e[2]) + log1p(-e[3]))
    e = BigFloat.(defects)
    all(z -> 0 < z < 1, e) || return 0.0
    bound = BigFloat(Inf)
    for j in 0:12
        z = zq(e)
        0 < z < 1//2 || return 0.0
        bound = min(bound, (log1p(-z)-log(z))/(big(2)^j*beta))
        if j == 12
            maximum(e) <= 1//1000 || return 0.0
            return Float64(min(bound, -log(25maximum(e))/(big(2)^j*beta)))
        end
        a = f(z)
        a < 1 || return 0.0
        b = f(zq((e[2], e[3], a)))
        b < 1 || return 0.0
        c = f(zq((e[3], a, b)))
        c < 1 || return 0.0
        e = (a, b, c)
    end
end

# The direct route uses the same spectra but evaluates its spatial input at 2β.
function approximate_direct_defects(spectra, beta)
    f(z) = z^2 / (2*(1-z)^2)
    a = approximate_logz(spectra["4x2"], 2beta)
    b = approximate_logz(spectra["4x4"], 2beta)
    c = approximate_logz(spectra["4x4"], beta)
    ex, et = -expm1(b-2a), -expm1(b-2c)
    x = f(ex)
    x < 1 || return (2.0, 2.0, 2.0)
    y = f(-expm1(2log1p(-ex) + log1p(-x)))
    y < 1 || return (2.0, 2.0, 2.0)
    t = f(-expm1(4log1p(-et) + 2log1p(-x) + log1p(-y)))
    return (x, y, t)
end

all(==("--direct"), ARGS) || error("usage: refine_scan.jl [--direct]")
direct = "--direct" in ARGS
records = Dict{String,Any}()
for g in (direct ? [p.g for p in POINTS] : [18, 20, 24])
    manifest = read_result(Certificate.manifest_path(root, g))
    spectra = Certificate.load_spectra(root, manifest, g)
    rows = []
    for k in (direct ? (40:320) : (120:240))
        beta = k/(40g)
        if direct
            e = approximate_direct_defects(spectra, beta)
        else
            a, b, c = [approximate_logz(spectra[s], beta) for s in ("2x2", "4x2", "4x4")]
            d = approximate_logz(spectra["4x4"], 2beta)
            e = (-expm1(b-2a), -expm1(c-2b), -expm1(d-2c))
        end
        push!(rows, Dict("beta_g" => k/40, "beta" => beta,
            "defects" => collect(e), "bound" => approximate_bound(e, direct ? 2beta : beta)))
    end
    best = rows[argmax([r["bound"] for r in rows])]
    records[string(g)] = Dict("best" => best, "rows" => rows)
    println("g=$g: $best")
end
name = direct ? "direct_temperature_scan.json" : "refinement_scan.json"
write_result(joinpath(root, "results", name),
    Dict("status" => "floating-point diagnostics, not certified", "fields" => records))
