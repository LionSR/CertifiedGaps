# Heuristic extrapolation (NOT certified) of the numerical scan: fits the ground-energy defect
# D(n) = A exp(-n/ell) to n = 3, 4 cells and an independent-excitation thermal prefactor, then
# predicts the first admissible even seed. Reads the gitignored results/numerical_delta_*.json
# written by scripts/numerical_scan.jl and writes the tracked results/numerical_summary.json.
using CertifiedGaps, CertifiedGaps.Results, CertifiedGaps.Seeds, JSON3, LinearAlgebra
BLAS.set_num_threads(1)
const MODEL_DIR = normpath(joinpath(@__DIR__, ".."))

"Summarize every stored numerical scan into results/numerical_summary.json."
function summarize()
    rows=Any[]
    paths=sort(filter(p->startswith(basename(p),"numerical_delta_")&&endswith(p,".json"),
        readdir(joinpath(MODEL_DIR,"results");join=true)))
    for path in paths
        data=JSON3.read(read(path,String)); diag=data.diagnostics
        n3=only(filter(r->r.n0==3,diag)); n4=only(filter(r->r.n0==4,diag))
        xi=1/log(n3.energy_difference/n4.energy_difference)
        amplitude=n4.energy_difference*exp(4/xi)
        gap=data.reference_gap
        probe=only(filter(r->r.n0==4&&r.beta0==8.0,data.scan))
        target=1-probe.eps_t
        # Fit a per-cell independent-excitation prefactor to the actual 16-spin purity.
        purity(a)=exp(8*(log1p(a*exp(-16gap))-2log1p(a*exp(-8gap))))
        lo=0.0; hi=100.0
        @assert purity(lo)>=target>=purity(hi)
        for _ in 1:80
            mid=(lo+hi)/2
            if purity(mid)>target; lo=mid; else; hi=mid; end
        end
        a=(lo+hi)/2
        model(n,b)=(-expm1(-b*amplitude*exp(-n/xi)),
            -expm1(2n*(log1p(a*exp(-2b*gap))-2log1p(a*exp(-b*gap)))))
        predicted=nothing
        for n in 2:2:40
            best=best_seed(scan_seeds(model,[n],0.25:0.25:64))
            if best!==nothing
                predicted=Dict("n0_cells"=>n,"beta0"=>best.β₀,"surrogate_bound"=>best.bound)
                break
            end
        end
        parts=split(data.delta,"//"); delta=parse(Float64,parts[1])/parse(Float64,parts[2])
        row=Dict("delta"=>data.delta,"reference_gap"=>gap,"reference_spins"=>data.reference_spins,
            "best_even_seed"=>data.best_even_seed,"ratio"=>data.bound_to_finite_ring_gap_ratio,
            "diagnostics"=>diag,"spectrum_seconds"=>data.seconds_spectra,
            "effective_energy_difference_decay_cells"=>xi,
            "decay_times_delta_two_thirds"=>xi*delta^(2/3),
            "energy_difference_amplitude"=>amplitude,"thermal_prefactor_per_cell"=>a,
            "surrogate_first_even_seed"=>predicted)
        push!(rows,row)
        println("delta=$(data.delta) gap=$gap best=$(data.best_even_seed) ",
            "ratio=$(data.bound_to_finite_ring_gap_ratio) seconds=$(data.seconds_spectra)")
        println("diagnostics=$(data.diagnostics)")
        println("xi=$xi scaled=$(xi*delta^(2/3)) prediction=$predicted")
    end
    root=normpath(joinpath(MODEL_DIR,"..",".."))
    write_result(joinpath(MODEL_DIR,"results","numerical_summary.json"),Dict(
        "metadata"=>run_metadata(Dict("fit"=>
                "two-point exponential ground-energy defect plus independent excitations");
            source_sha256=source_hashes(root;dirs=(),
                files=("models/dimerized_heisenberg/scripts/summarize_numerics.jl",)),
            status="heuristic_extrapolation_not_certified",certified=false),
        "input_hashes"=>Dict(relpath(p,root)=>sha256_file(p) for p in paths),
        "rows"=>rows,"limitations"=>"Finite-ring gaps, short-size exponential fits, and " *
            "independent excitation ansatz; extrapolated seeds are not validated or " *
            "guaranteed. Odd n=3 used only for fitting, never as theorem seed."))
end
summarize()
