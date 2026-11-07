# Complete floating-point spectra and seed diagnostics at moderate transverse fields.
using LinearAlgebra
using CertifiedGaps.Results
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
include(joinpath(@__DIR__, "..", "src", "SquareIsing.jl"))
using .SquareIsing

BLAS.set_num_threads(1)
root = dirname(@__DIR__)
torus_sizes = ((2, 2), (4, 2), (4, 4))
fields = isempty(ARGS) ? (4, 5, 6, 8) : parse.(Int, ARGS)

function logz(energies, beta)
    e0 = first(energies)
    return -beta*e0 + log(sum(exp(-beta*(e-e0)) for e in energies))
end

function recurrence_diagnostic(defects)
    f(z) = z^2 / (2*(1-z)^2)
    e = defects
    for j in 0:24
        all(z -> 0 <= z < 1, e) || return false
        maximum(e) <= 0.001 && return true
        a = f(1-(1-e[1])^4*(1-e[2])^2*(1-e[3]))
        a < 1 || return false
        b = f(1-(1-e[2])^4*(1-e[3])^2*(1-a))
        b < 1 || return false
        c = f(1-(1-e[3])^4*(1-a)^2*(1-b))
        e = (a, b, c)
    end
    return false
end

for g in fields
    spectra = []
    finite = Dict{String,Any}()
    for (n, m) in torus_sizes
        seconds = @elapsed s = float_ring_spectrum(row_model(m, g), n)
        push!(spectra, s.E)
        finite["$(n)x$m"] = Dict("sites" => string(n*m), "ground_energy" => s.E[1],
            "energy_per_site" => s.E[1]/(n*m), "gap" => s.E[2]-s.E[1])
        println("g=$g $(n)x$m: $(finite["$(n)x$m"]); $seconds seconds")
        flush(stdout)
    end
    rows = []
    for k in 20:640
        beta = k/(40g)
        a, b, c = [logz(s, beta) for s in spectra]
        d = logz(spectra[3], 2beta)
        e = (-expm1(b-2a), -expm1(c-2b), -expm1(d-2c))
        q = exp(4*(b-2a)+2*(c-2b)+(d-2c))
        push!(rows, Dict("beta" => beta, "beta_g" => k/40, "defects" => collect(e),
            "cube_lower" => q, "recurrence_closes" => recurrence_diagnostic(e)))
    end
    best = rows[argmax([r["cube_lower"] for r in rows])]
    result = Dict("status" => "floating-point diagnostics; no uniform gap certificate",
        "g" => string(g), "finite_volume" => finite, "seed_scan" => rows,
        "largest_cube_lower" => best,
        "closing_seeds" => string(count(r -> r["recurrence_closes"], rows)),
        "ground_energy_mismatches" => [spectra[2][1]-2spectra[1][1],
            spectra[3][1]-2spectra[2][1]])
    write_result(joinpath(root, "results", "moderate_g$g.json"), result)
    println("g=$g: best Q=$best; closing seeds=$(result["closing_seeds"])")
    flush(stdout)
end
