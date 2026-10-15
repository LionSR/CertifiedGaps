using CertifiedGaps.Results, JSON3, Printf

"Preserve every fresh scan row and render the requested ratio curve with PGFPlots."
function report_transfer()
    root = @__DIR__
    results = joinpath(root, "results")
    cert = parse_exact_json(read(joinpath(root, "certificate.json"), String))
    cases = cert["cases"]
    scans, curve = Any[], Any[]
    Ds = sort(unique(vcat(collect(1:1//4:4), [641//200, 1603//500])))
    for D in Ds
        stem = "D$(numerator(D))_$(denominator(D))"
        path = joinpath(results, "summary_$stem.json")
        data = JSON3.read(read(path, String), Dict{String,Any})
        finepath = joinpath(results, "beta_refinement_$stem.json")
        fine = isfile(finepath) ? JSON3.read(read(finepath, String), Dict{String,Any}) : nothing
        best = fine === nothing ? data["best"] : fine["best"]
        push!(scans, Dict("D" => string(D), "coarse" => data, "fine" => fine))
        gap = only(filter(x -> x["L"] == 12, data["sizes"]))["gap"]
        push!(curve, Dict("D" => Float64(D), "D_exact" => string(D),
            "bound" => best === nothing ? nothing : best["bound"], "gap_L12" => gap,
            "ratio" => best === nothing ? nothing : best["bound"]/gap,
            "status" => "numerical screening, n0<=6; missing means no tested seed"))
    end
    sources = source_hashes(root; dirs = (),
        files = ("model.jl", "scan.jl", "refine_beta.jl", "report_transfer.jl"))
    write_result(joinpath(results, "fresh_scan.json"), Dict("certified" => false,
        "scans" => scans, "source_hashes" => sources,
        "run_metadata" => run_metadata(Dict("grid" => "1:1/4:4; 3.205,3.206 refined");
            source_sha256 = sources, status = "floating-point screening only")))
    chosen = [only(filter(c -> c["id"] == id, cases)) for id in
        ("transport-threshold-D1603//500", "space-lift-1-D4//1")]
    exactvalue(x) = Float64(parse(Rational{BigInt}, x))
    for c in chosen
        push!(curve, Dict("D" => exactvalue(c["D_interval"][1]),
            "D_exact" => c["D_interval"][1], "bound" => exactvalue(c["bound"]),
            "gap_L12_interval" => c["finite_gap_L12"],
            "ratio_interval" => c["ratio_to_gap_L12"],
            "status" => "certified numerator and finite L=12 denominator interval"))
    end
    write_result(joinpath(results, "transfer_curve.json"), Dict("rows" => curve,
        "denominator" => "finite-ring Delta_12; no thermodynamic gap asserted",
        "source_hashes" => sources,
        "certificate_sha256" => sha256_file(joinpath(root, "certificate.json"))))
    open(joinpath(results, "transfer_curve.csv"), "w") do io
        println(io, "D,bound,ratio_lower,ratio_upper,status")
        for r in curve
            if haskey(r, "ratio_interval")
                l, h = exactvalue.(r["ratio_interval"])
            else
                l = h = r["ratio"] === nothing ? "" : r["ratio"]
            end
            println(io, join((r["D"], something(r["bound"], ""), l, h, r["status"]), ','))
        end
    end
    numerical = filter(r -> haskey(r, "ratio") && r["ratio"] !== nothing, curve)
    points = join(("($(r["D"]),$(r["ratio"]))" for r in numerical), "\n")
    certified = join(("($(exactvalue(c["D_interval"][1]))," *
        "$(exactvalue(c["ratio_to_gap_L12"][1])))" for c in chosen), "\n")
    tex = raw"""
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
""" * certified * raw"""
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
    path = joinpath(results, "transfer_curve.tex")
    write(path, tex)
    command = `pdflatex -interaction=nonstopmode -halt-on-error -output-directory=$results $path`
    run(pipeline(command;
        stdout = joinpath(root, "plot.log"), stderr = stderr))
    run(`pdftocairo -svg $(joinpath(results, "transfer_curve.pdf"))
        $(joinpath(results, "transfer_curve.svg"))`)
    run(`pdftocairo -png -singlefile -scale-to 1500 $(joinpath(results, "transfer_curve.pdf"))
        $(joinpath(results, "transfer_curve"))`)
end

if abspath(PROGRAM_FILE) == @__FILE__
    report_transfer()
end
