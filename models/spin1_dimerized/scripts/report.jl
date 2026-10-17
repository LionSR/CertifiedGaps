# Rebuilds results/summary.json and the "## Results from completed runs" section of METHOD.md
# from the scans, the certificates and the regenerated spectra (scripts/scan.jl, all c6 files).
include(joinpath(dirname(@__DIR__), "src", "DimerizedCertificate.jl"))
using .DimerizedCertificate, CertifiedGaps.Results, CertifiedGaps.Seeds, JSON3, Printf
const D = DimerizedCertificate
const DELTAS, ROOT, ratstr, charge_groups = D.DELTAS, D.ROOT, D.str, D.charge_groups
"Certificate summary: exact chain bound, interval bound, seed and run time."
function certificate_entry(cert, gap)
    bound = parse(Float64, string(cert.gap_lower_bound_decimal))
    Dict("n0" => cert.seed.n0, "beta0" => cert.seed.beta0, "bound" => bound,
        "interval_bound" => cert.interval_bound_decimal, "seconds" => cert.seconds,
        "ratio_to_12_spin_gap_numerical" => bound / gap)
end
lower12(x) = let r=Rational{BigInt}(x), d=big(10)^12
    q=fld(numerator(r)*d,denominator(r))
    string(q÷d,".",lpad(string(q%d),12,Char(48)))
end
"Numerical two-point finite-size extrapolation, never an input to certification."
function estimate_next(scan,gap)
    scan.best!==nothing && return nothing
    s=sort(collect(scan.summary);by=x->x.n0)
    length(s)==2 || return nothing
    a,b=s[1].energy_defect,s[2].energy_defect
    (a>b>0 && gap>0) || return nothing
    rate=log(a/b); amplitude=a*exp(2rate)
    for n in 4:200
        beta=ceil(4log(12n/0.05)/gap)/4
        es=1-exp(-beta*amplitude*exp(-rate*n));et=12n*exp(-beta*gap)
        bound=certified_gap_bound(es,et,beta).bound
        if es<=0.05 && bound!==nothing
            return Dict("n0"=>n,"beta0"=>beta,"decay_rate"=>rate,"amplitude"=>amplitude,"assumed_gap"=>gap,"predicted_epsilon_s"=>es,"predicted_epsilon_t"=>et,"predicted_bound_numerical_only"=>bound)
        end
    end
    nothing
end
rows=Any[]
for delta in DELTAS
    tag=Int(10delta)
    scan=JSON3.read(read(joinpath(ROOT,"results","scan_d$tag.json"),String))
    spec=JSON3.read(read(joinpath(ROOT,"results","spectrum_d$(tag)_c6.json"),String))
    cp=joinpath(ROOT,"results","certificate_d$tag.json")
    cert=isfile(cp) ? JSON3.read(read(cp,String)) : nothing
    push!(rows,Dict("delta"=>ratstr(delta),"best_numerical_seed"=>scan.best,"purity_and_energy_scaling"=>scan.summary,"gap_12_spins_numerical"=>spec.gap,"certificate"=>cert===nothing ? nothing : certificate_entry(cert,spec.gap),"extrapolation_numerical_only"=>estimate_next(scan,spec.gap)))
end
metadata=model_metadata(ROOT,
    Dict("purpose"=>"parameter comparison and numerical extrapolation");D.SOURCES...)
write_result(joinpath(ROOT,"results","summary.json"),Dict("metadata"=>metadata,"parameters"=>rows))
reportpath=joinpath(ROOT,"METHOD.md")
base=first(split(read(reportpath,String),"\n## Results from completed runs"))
open(reportpath,"w") do io
    write(io,base)
    println(io,"\n## Results from completed runs\n")
    println(io,"Best numerical seed per parameter, over n=2,3 and beta=0.25:0.25:32:\n")
    println(io,"| delta | n | beta | epsilon_s | epsilon_t | numerical recurrence bound | ED gap (12 spins) |")
    println(io,"|---|---:|---:|---:|---:|---:|---:|")
    for row in rows
        b=row["best_numerical_seed"]
        if b===nothing
            @printf(io,"| %s | — | — | — | — | no seed | %.9f |\n",row["delta"],row["gap_12_spins_numerical"])
        else
            @printf(io,"| %s | %d | %.2f | %.9f | %.9f | %.9f | %.9f |\n",row["delta"],b.n0,b.beta0,b.epsilon_s,b.epsilon_t,b.bound,row["gap_12_spins_numerical"])
        end
    end
    println(io,"\nVerified certificates (the 0.9, n=2 seed was chosen before the larger-ring scan):\n")
    println(io,"| delta | n | exact beta | certified lower bound | ratio to 12-spin ED gap | verification seconds | validity |")
    println(io,"|---|---:|---:|---:|---:|---:|---|")
    for row in rows
        c=row["certificate"];c===nothing && continue
        @printf(io,"| %s | %d | %s | %s | %.6f | %.2f | even L >= %d cells; N divisible by 4, N >= %d spins |\n",row["delta"],c["n0"],c["beta0"],lower12(c["bound"]),c["ratio_to_12_spin_gap_numerical"],c["seconds"],c["n0"],4cld(c["n0"],2))
    end
    println(io,"\nThe table gives each certificate's exact rational chain bound `gap_lower_bound`; its interval bound `certified_gap_bound` is `interval_bound_decimal` in the certificate (and in README.md). Table bounds are rounded down to twelve decimal places.\n")
    println(io,"| delta | max min purity, n=2 | max min purity, n=3 | E0(4)-2E0(2) | E0(6)-2E0(3) | extrapolated (n,beta), numerical only |")
    println(io,"|---|---:|---:|---:|---:|---|")
    for row in rows
        s=sort(collect(row["purity_and_energy_scaling"]);by=x->x.n0)
        e=row["extrapolation_numerical_only"]
        ex=e===nothing ? "not needed" : "($(e["n0"]), $(e["beta0"]))"
        @printf(io,"| %s | %.9f | %.9f | %.10f | %.10f | %s |\n",row["delta"],s[1].max_min_purity,s[2].max_min_purity,s[1].energy_defect,s[2].energy_defect,ex)
    end
    println(io,"\nThe extrapolation fits delta(n)=A exp(-kappa n) through n=2,3 and holds the numerical twelve-spin gap fixed. It estimates epsilon_t=12n exp(-beta gap) and epsilon_s=1-exp(-beta delta(n)), targeting both below 0.05. Two sizes cannot validate this scaling or an infinite-volume gap; these are planning estimates, not demonstrated seeds.\n")
    println(io,"| cells | spins | full dimension | largest numerical block | largest full Sz block |")
    println(io,"|---:|---:|---:|---:|---:|")
    for cells in (2,3,4,6)
        data=JSON3.read(read(joinpath(ROOT,"results","spectrum_d8_c$cells.json"),String))
        @printf(io,"| %d | %d | %d | %d | %d |\n",cells,2cells,9^cells,maximum(s.dimension for s in data.sectors),length(charge_groups(cells)[1]))
    end
    println(io,"\nFull sector dimension lists and multiplicities are in `results/certificate_d*.json` ",
        "and in the regenerated `results/spectrum_d*_c*.json` (`scripts/scan.jl`).\n")
    println(io,"| delta | dense search seconds (sum across sizes) | twelve-spin seconds |")
    println(io,"|---|---:|---:|")
    for delta in DELTAS
        specs=[JSON3.read(read(joinpath(ROOT,"results","spectrum_d$(Int(10delta))_c$c.json"),String)) for c in (2,3,4,6)]
        @printf(io,"| %s | %.2f | %.2f |\n",ratstr(delta),sum(s.seconds for s in specs),specs[end].seconds)
    end
end
