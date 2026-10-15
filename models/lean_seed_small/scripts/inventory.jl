# NON-CERTIFYING inventory of existing spectral data and the exact Lean seed targets.
using JSON, SHA, Printf
const ROOT = joinpath(@__DIR__, "..", "results")
const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
rat(s::AbstractString) = occursin("//", s) ?
    (a = split(s, "//"); parse(BigInt, a[1]) // parse(BigInt, a[2])) :
    rationalize(BigInt, parse(BigFloat, s); tol=0)
num(x::AbstractString) = Float64(rat(x))
num(x::Number) = Float64(x)
function readjson(rel)
    p = joinpath(REPO, "models", rel)
    JSON.parse(endswith(p, ".gz") ? read(`gzip -dc $p`, String) : read(p, String))
end
function levels(blocks; key="eigenvalues")
    reduce(vcat, [repeat(num.(b[key]), get(b, "multiplicity", 1)) for b in blocks])
end
function readcompact(rel)
    open(joinpath(REPO, "models", rel)) do io
        @assert read(io, 8) == b"SPIN1LD1"
        nr = Int(ltoh(read(io, UInt64))); nt = Int(ltoh(read(io, UInt64)))
        read!(io, Vector{Float64}(undef, nt))
        es = ltoh.(read!(io, Vector{Float64}(undef, nr)))
        read!(io, Vector{UInt8}(undef, nr))
        counts = read!(io, Vector{UInt8}(undef, nr))
        @assert eof(io)
        reduce(vcat, [fill(es[i], Int(counts[i])) for i in eachindex(es)])
    end
end
logpart(es, beta) = -beta * minimum(es) + log(sum(exp.(-beta .* (es .- minimum(es)))))
rows = []
function addrow(name, n, q, beta, es, et, source_s, source_t, E, path)
    @assert length(E) == q^(2n)
    # The Lean target minus the stored external enclosure is a conservative planning slack,
    # conditional on that enclosure. It is NOT assumed to be true in any Lean proof.
    slack = min(es - num(source_s), et - num(source_t))
    @assert slack > 0
    eta = min(1e-6, slack / 16)
    e0 = minimum(E)
    cutoff = (log(length(E)) + log(1 / eta)) / beta
    k = count(e -> e - e0 <= cutoff, E)
    push!(rows, Dict("model"=>name, "n0"=>n, "local_dimension"=>q,
        "dimension"=>length(E), "beta"=>beta, "eps_s"=>es, "eps_t"=>et,
        "enclosure_slack"=>slack, "tail_budget"=>eta,
        "uniform_eigen_radius_budget"=>eta/(4beta), "cutoff_above_ground"=>cutoff,
        "levels_below_cutoff"=>k, "ground_energy"=>e0,
        "actual_weight_above_cutoff"=>sum(exp(-beta*(e-e0)) for e in E if e-e0>cutoff),
        "source"=>path, "source_hash"=>bytes2hex(sha256(read(joinpath(REPO,"models",path)))),
        "status"=>"NON-CERTIFYING Float64 planning estimate"))
end
for (tag, b, es, et) in [("4_1",2.25,.06411,.04053),
                         ("1603_500",3.048,.19164415,.046462453)]
    c=readjson("spin1_largeD/results/certificate_D$tag.json")
    p="spin1_largeD/results/enclosed_D$(tag)_L12.levels.bin"
    addrow("LargeD $tag",6,3,b,es,et,c["eps_s"],c["eps_t"],readcompact(p),p)
end
let p="heisenberg_ladder/results/verified_L8.json.gz"
    c=readjson("heisenberg_ladder/results/certificate.json"); d=readjson(p)
    # This artifact uses the same names as the original model certificate.
    addrow("Heisenberg ladder J4",4,4,2.,.03509,.02463,
        c["eps_s"],c["eps_t"],num.(d["energies"]),p)
end
for (tag, b, es, et) in [("1_1",5.25,.0006937,.0007327),
    ("3_4",5.25,.002343,.003123),("1_2",5.5,.01075,.01116),
    ("7_20",5.75,.03266,.02744),("1_4",6.25,.07987,.04477)]
    p="heisenberg_staggered_field/results/certificate_$tag.json"; c=readjson(p)
    addrow("Staggered h$tag",8,2,b,es,et,c["epsilon_s"],c["epsilon_t"],
        levels(c["large_sectors"]),p)
end
for (tag,b,es,et) in [("1_2",6.,.01133,.01224),("3_10",6.5,.09516,.0518),
                      ("1_4",7.25,.1658,.05123)]
    p="dimerized_heisenberg/results/verified_delta_$tag.json.gz"; c=readjson(p)
    addrow("Dimerized half d$tag",4,4,b,es,et,c["epsilon_s_upper"],
        c["epsilon_t_upper"],levels(c["rings"][2]["blocks"]),p)
end
let p="spin1_dimerized/results/spectrum_d9_c4.json"
    c=readjson("spin1_dimerized/results/certificate_d9.json")
    addrow("Dimerized one d9/10",2,9,3.5,.04972,.03206,c["epsilon_s"],
        c["epsilon_t"],num.(readjson(p)["energies"]),p)
end
let p="two_leg_tfim/results/certificate_g3_1.json.gz"
    c=readjson(p)
    addrow("Ising ladder G3",4,4,1.5,.1125,.05162,c["epsilon_s"],
        c["epsilon_t"],num.(c["spectra"][2]["energies"]),p)
end
# Pull each rounded seed bound directly from the current Lean source, not from README tables.
src=read(joinpath(REPO,"lean/CertifiedGaps/Models/IsingLongitudinal.lean"),String)
pat=r"⟨([^,]+), ([^,]+), ([^,]+), \{\s*s := \[([^,]+),[\s\S]*?t := \[([^,]+),"
seeds=readjson("tfim_longitudinal/results/seeds.json")["points"]
ms=collect(eachmatch(pat,src)); @assert length(ms)==length(seeds)==30
for (m, s) in zip(ms,seeds)
    for (capture,key) in [(m[1],"g"),(m[2],"h"),(m[3],"beta0")]
        a=split(replace(strip(capture)," "=>""),"/")
        q=length(a)==1 ? parse(BigInt,a[1])//big(1) : parse(BigInt,a[1])//parse(BigInt,a[2])
        @assert q==rat(s[key])
    end
    p="tfim_longitudinal/results/"*s["file"]; c=readjson(p)["certificate"]
    es=parse(Float64,strip(m[4])); et=parse(Float64,strip(m[5]))
    addrow("Ising g$(s["g"]) h$(s["h"])",6,2,num(s["beta0"]),es,et,
        s["epsilon_s"],s["epsilon_t"],levels(c["large"];key="values"),p)
end
open(joinpath(ROOT,"inventory.json"),"w") do io
    JSON.print(io,Dict("status"=>"NON-CERTIFYING", "rows"=>rows))
end
open(joinpath(ROOT,"inventory.tsv"),"w") do io
    println(io,"model\tdimension\tbeta\teps_s\teps_t\tslack\ttail_budget\teigen_radius\tcutoff\tk")
    for r in rows
        @printf(io,"%s\t%d\t%.6g\t%.9g\t%.9g\t%.3g\t%.3g\t%.3g\t%.5g\t%d\n",
            r["model"],r["dimension"],r["beta"],r["eps_s"],r["eps_t"],
            r["enclosure_slack"],r["tail_budget"],r["uniform_eigen_radius_budget"],
            r["cutoff_above_ground"],r["levels_below_cutoff"])
    end
end
println("NON-CERTIFYING inventory: $(length(rows)) parameter points")
