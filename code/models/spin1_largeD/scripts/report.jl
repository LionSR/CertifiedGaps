# Summarise the screening scan and draw the bound/gap curve (floating point, not a proof).
#
#   julia --project=julia models/spin1_largeD/scripts/report.jl [--certified-only]
#
# Inputs: the gitignored results/summary_D*.json and results/beta_refinement_D*.json from
# scripts/scan.jl and scripts/refine_beta.jl, and results/certificate.json. Outputs:
#   results/scan_summary.json  per-D sizes, purity diagnostics and best seeds (no raw rows);
#   results/transfer_curve.json  numerical ratios and certified exact ratio intervals;
#   results/transfer_curve.pdf  rendered with pdflatex/PGFPlots in a temporary directory.
using CertifiedGaps.Results, JSON3

const ROOT = dirname(@__DIR__)
const RESULTS = joinpath(ROOT, "results")
const SCANNED = sort(unique(vcat(collect(1:1//4:4), [641//200, 1603//500])))
const CERTIFIED_IDS = ("transport-D1603//500", "space-lift-1-D4//1")

readjson(path) = JSON3.read(read(path, String), Dict{String,Any})
exactvalue(x) = Float64(parse(Rational{BigInt}, x))

"Per-D scan summary: coarse best, refined best and the purity/delta diagnostics."
function scan_summary(D)
    stem = "D$(numerator(D))_$(denominator(D))"
    coarse = readjson(joinpath(RESULTS, "summary_$stem.json"))
    finepath = joinpath(RESULTS, "beta_refinement_$stem.json")
    fine = isfile(finepath) ? readjson(finepath) : nothing
    return Dict("D" => string(D), "sizes" => coarse["sizes"],
        "diagnostics" => coarse["diagnostics"], "best_coarse" => coarse["best"],
        "best_fine" => fine === nothing ? nothing : fine["best"],
        "seconds" => sum(x["seconds"] for x in coarse["sizes"]))
end

"PGFPlots source for the numerical curve and the certified markers."
function curve_tex(numerical, certified)
    points = join(("($(r["D"]),$(r["ratio"]))" for r in numerical), "\n")
    marks = join(("($(r["D"]),$(exactvalue(r["ratio_interval"][1])))" for r in certified),
        "\n")
    return raw"""
\documentclass[border=10pt]{standalone}
\usepackage{pgfplots}
\pgfplotsset{compat=1.18}
\begin{document}
\begin{tikzpicture}
\begin{axis}[width=15cm,height=8cm,xmin=1,xmax=4.08,ymin=0,ymax=.43,
  xlabel={$D$},ylabel={Uniform bound $B\,/\,\Delta_{12}$},
  title={Spin-1 large-$D$ chain: certified bounds and ED screening},
  xtick={1,1.5,2,2.5,3,3.5,4},ytick={0,.1,.2,.3,.4},
  grid=major,grid style={gray!18},axis lines=left,clip=false,
  legend style={at={(.02,.98)},anchor=north west,draw=none,font=\small},
  tick label style={font=\small},label style={font=\small}]
\addplot[blue!65!black,thick,mark=o,mark options={fill=white}] coordinates {
""" * points * raw"""
};
\addlegendentry{Numerical scan, $n_0\le6$}
\addplot[orange!80!black,only marks,mark=diamond*,mark size=3.2pt] coordinates {
""" * marks * raw"""
};
\addlegendentry{Certified bounds / enclosed $\Delta_{12}$}
\node[anchor=west,font=\small,align=left,text=gray!75!black]
  at (axis cs:1.2,.18)
  {Sampled $D\le3.205$: no seed on the tested grid.\\This does not imply a zero gap.};
\node[anchor=north east,font=\small,align=right] at (axis cs:3.98,.388)
  {$D=4$: transfer seed $n_0=12$\\plus finite bridge to $L\ge6$};
\node[anchor=north,font=\small,align=center] at (rel axis cs:.5,-.17)
  {Denominator: finite periodic $L=12$ gap, not the infinite-chain gap.\\
   Curves connect sampled numerical values only; certificates are the marked points.};
\end{axis}
\end{tikzpicture}
\end{document}
"""
end

"The certified points of the curve, copied from the cases of results/certificate.json."
function certified_points()
    cases = read_exact_json(joinpath(RESULTS, "certificate.json"))["cases"]
    return [let c = only(filter(c -> c["id"] == id, cases))
                Dict("D" => exactvalue(c["D_interval"][1]), "id" => id,
                    "D_interval" => c["D_interval"], "bound" => c["bound"],
                    "gap_L12_interval" => c["finite_gap_L12"],
                    "ratio_interval" => c["ratio_to_gap_L12"],
                    "status" => "certified numerator and finite L=12 denominator interval")
            end for id in CERTIFIED_IDS]
end

"""
    report(; certified_only = false)

Write the scan summary and the curve data, then render the PDF. With `certified_only`, keep the
stored numerical scan (results/scan_summary.json and the curve's numerical points) and refresh
only the certified points from results/certificate.json.
"""
function report(; certified_only = false)
    metadata = model_metadata(ROOT, Dict("scanned_D" => string.(SCANNED),
            "certified_only" => certified_only); dirs = (),
        files = ("src/model.jl", "scripts/scan.jl", "scripts/refine_beta.jl",
            "scripts/report.jl"))
    if certified_only
        numerical = readjson(joinpath(RESULTS, "transfer_curve.json"))["numerical"]
    else
        scans = [scan_summary(D) for D in SCANNED]
        write_result(joinpath(RESULTS, "scan_summary.json"), Dict("certified" => false,
            "grid" => "D = 1:1/4:4 and 3.205, 3.206; n0 = 2..6; beta = 0.125:0.125:20; " *
                      "refined beta = 2.75:0.002:3.75 at n0 = 6 for D = 3.205, 3.206",
            "note" => "floating-point screening; odd n0 rows are diagnostics only",
            "scans" => scans, "metadata" => metadata))
        numerical = Any[]
        for s in scans
            best = something(s["best_fine"], s["best_coarse"], Some(nothing))
            gap = only(filter(x -> x["L"] == 12, s["sizes"]))["gap"]
            push!(numerical, Dict("D" => exactvalue(s["D"]), "D_exact" => s["D"],
                "bound" => best === nothing ? nothing : best["bound"], "gap_L12" => gap,
                "ratio" => best === nothing ? nothing : best["bound"] / gap,
                "status" => "numerical screening, n0<=6; missing means no tested seed"))
        end
    end
    certified = certified_points()
    write_result(joinpath(RESULTS, "transfer_curve.json"), Dict(
        "numerical" => numerical, "certified" => certified,
        "denominator" => "finite-ring Delta_12; no thermodynamic gap asserted",
        "metadata" => metadata,
        "certificate_sha256" => sha256_file(joinpath(RESULTS, "certificate.json"))))
    shown = filter(r -> r["ratio"] !== nothing, numerical)
    mktempdir() do dir
        tex = joinpath(dir, "transfer_curve.tex")
        write(tex, curve_tex(shown, certified))
        run(pipeline(`pdflatex -interaction=nonstopmode -halt-on-error
            -output-directory=$dir $tex`; stdout = devnull))
        cp(joinpath(dir, "transfer_curve.pdf"), joinpath(RESULTS, "transfer_curve.pdf");
            force = true)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    report(; certified_only = "--certified-only" in ARGS)
end
