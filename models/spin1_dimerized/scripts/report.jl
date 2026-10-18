# Rebuilds results/summary.json and the "## Results from completed runs" section of METHOD.md
# from the scans, the certificates and the cached scan spectra (run scripts/scan.jl for every δ
# of DELTAS first). Floating-point summaries except the certificate rows.
include(joinpath(dirname(@__DIR__), "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate, CertifiedGaps, CertifiedGaps.Results, JSON3, Printf
using CertifiedGaps.GapCertificates: certificate_row
using CertifiedGaps.ExactArithmetic: parse_rational
const D = DimerizedCertificate

readjson(path) = JSON3.read(read(path, String))
results(name) = joinpath(D.ROOT, "results", name)
spectrum(δ, cells) = readjson(results("spectrum_d$(Int(10δ))_c$cells.json"))

"Exact decimal `x` (a `num//den` string) rounded down to twelve places."
function lower12(x)
    r = parse_rational(x)
    q = fld(numerator(r) * big(10)^12, denominator(r))
    string(q ÷ big(10)^12, ".", lpad(string(q % big(10)^12), 12, '0'))
end

"""
Numerical two-point finite-size extrapolation, never an input to certification: fits
δ(n) = A exp(−κ n) through the energy mismatches at n = 2, 3, holds the twelve-spin gap fixed,
and returns the smallest n with predicted εₛ, εₜ ≤ 0.05 and an admissible bound.
"""
function estimate_next(scan, gap)
    scan.best === nothing || return nothing
    s = sort(collect(scan.summary); by = x -> x.n0)
    length(s) == 2 || return nothing
    a, b = s[1].energy_defect, s[2].energy_defect
    a > b > 0 && gap > 0 || return nothing
    rate = log(a / b)
    amplitude = a * exp(2rate)
    for n in 4:200
        β = ceil(4log(12n / 0.05) / gap) / 4
        εs = 1 - exp(-β * amplitude * exp(-rate * n))
        εt = 12n * exp(-β * gap)
        bound = certified_gap_bound(εs, εt, β).bound
        εs <= 0.05 && bound !== nothing && return Dict("n0" => n, "beta0" => β,
            "decay_rate" => rate, "amplitude" => amplitude, "assumed_gap" => gap,
            "predicted_epsilon_s" => εs, "predicted_epsilon_t" => εt,
            "predicted_bound_numerical_only" => bound)
    end
    return nothing
end

rows = map(D.DELTAS) do δ
    tag = Int(10δ)
    scan = readjson(results("scan_d$tag.json"))
    gap = spectrum(δ, 6).gap
    path = D.certificate_path(δ)
    cert = isfile(path) ? read_result(path) : nothing
    entry = cert === nothing ? nothing : merge(certificate_row(cert), Dict(
        "seconds" => cert["metadata"]["config"]["seconds"],
        "ratio_to_12_spin_gap_numerical" =>
            parse(Float64, string(cert["gap_lower_bound_decimal"])) / gap,
        "spectra" => cert["spectra"]))
    Dict("delta" => D.str(δ), "best_numerical_seed" => scan.best,
        "purity_and_energy_scaling" => scan.summary, "gap_12_spins_numerical" => gap,
        "certificate" => entry, "extrapolation_numerical_only" => estimate_next(scan, gap))
end
metadata = model_metadata(D.ROOT, Dict("purpose" => "parameter comparison and extrapolation"))
write_result(results("summary.json"), Dict("metadata" => metadata, "parameters" => rows))

reportpath = joinpath(D.ROOT, "METHOD.md")
base = first(split(read(reportpath, String), "\n## Results from completed runs"))
open(reportpath, "w") do io
    write(io, base)
    println(io, "\n## Results from completed runs\n")
    println(io, "Best numerical seed per parameter, over n = 2, 3 and β = 0.25:0.25:32 ",
        "(floating point):\n")
    println(io, "| delta | n | beta | epsilon_s | epsilon_t | numerical recurrence bound | ",
        "ED gap (12 spins) |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|")
    for row in rows
        b, g = row["best_numerical_seed"], row["gap_12_spins_numerical"]
        if b === nothing
            @printf(io, "| %s | — | — | — | — | no seed | %.9f |\n", row["delta"], g)
        else
            @printf(io, "| %s | %d | %.2f | %.9f | %.9f | %.9f | %.9f |\n", row["delta"],
                b.n0, b.beta0, b.eps_s, b.eps_t, b.bound, g)
        end
    end
    println(io, "\nVerified certificates (the 0.9, n=2 seed was chosen before the larger-ring ",
        "scan):\n")
    println(io, "| delta | n | exact beta | interval bound | exact chain bound | ",
        "ratio to 12-spin ED gap | certify seconds | validity |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|---|")
    for row in rows
        c = row["certificate"]
        c === nothing && continue
        n0 = c["n0"]
        @printf(io, "| %s | %s | %s | %.9f | %s | %.6f | %.1f | even L >= %s cells |\n",
            row["delta"], n0, c["beta0"], floor(c["interval_bound_decimal"] * 1e9) / 1e9,
            lower12(c["gap_lower_bound"]), c["ratio_to_12_spin_gap_numerical"],
            c["seconds"], n0)
    end
    println(io, "\nThe chain bounds are rounded down to twelve decimal places, the interval ",
        "bounds to nine.\n")
    println(io, "| delta | ring (cells) | dimension | blocks | largest block | max radius | ",
        "block seconds |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|")
    for row in rows
        c = row["certificate"]
        c === nothing && continue
        for cells in sort(parse.(Int, collect(keys(c["spectra"]))))
            s = c["spectra"][string(cells)]
            @printf(io, "| %s | %d | %d | %d | %d | %.2e | %.1f |\n", row["delta"], cells,
                s["dimension"], s["blocks"], s["largest_block"], s["max_radius_decimal"],
                s["seconds"])
        end
    end
    println(io, "\n| delta | max min purity, n=2 | max min purity, n=3 | E0(4)-2E0(2) | ",
        "E0(6)-2E0(3) | extrapolated (n,beta), numerical only |")
    println(io, "|---|---:|---:|---:|---:|---|")
    for row in rows
        s = sort(collect(row["purity_and_energy_scaling"]); by = x -> x.n0)
        e = row["extrapolation_numerical_only"]
        ex = e === nothing ? "not needed" : "($(e["n0"]), $(e["beta0"]))"
        @printf(io, "| %s | %.9f | %.9f | %.10f | %.10f | %s |\n", row["delta"],
            s[1].max_min_purity, s[2].max_min_purity, s[1].energy_defect, s[2].energy_defect, ex)
    end
    println(io, "\nThe extrapolation fits delta(n)=A exp(-kappa n) through n=2,3 and holds the ",
        "numerical twelve-spin gap fixed. It estimates epsilon_t=12n exp(-beta gap) and ",
        "epsilon_s=1-exp(-beta delta(n)), targeting both below 0.05. Two sizes cannot validate ",
        "this scaling or an infinite-volume gap; these are planning estimates, not demonstrated ",
        "seeds.\n")
    println(io, "| delta | scan seconds (sum across sizes) | twelve-spin seconds |")
    println(io, "|---|---:|---:|")
    for δ in D.DELTAS
        t = [spectrum(δ, c).seconds for c in (2, 3, 4, 6)]
        @printf(io, "| %s | %.2f | %.2f |\n", D.str(δ), sum(t), t[end])
    end
end
