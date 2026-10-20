# Twenty-spin numerical diagnostics (NOT certified). Every Sz, cell-momentum and spin-flip sector
# of the 20-spin ring is enumerated; large blocks keep 40 Krylov Ritz values, so the spectrum is
# partial and nothing here enters a certificate. The complete spectra of the smaller rings are
# the core `float_ring_spectrum` of `rotated_model` (same spectrum as the physical ring).
# `run_diagnostic` writes the gitignored results/diagnostic20_<p>_<q>.json (about 390 kB, mostly
# Ritz values); `summarize_diagnostics` condenses the six files into the tracked
# results/diagnostic20_summary.json.
module Diagnostics
using CertifiedGaps, CertifiedGaps.Results, KrylovKit, LinearAlgebra, Random, SparseArrays
using CertifiedGaps.ExactDiagonalization: PartialSpectrum, log_partition, logZ_bounds
using CertifiedGaps.VerifiedSpectra: float_ring_spectrum
using ..DimerizedHeisenberg: rotated_model, MODEL_DIR

"Translate a configuration of `L` spins by one two-spin cell."
rotate_cell(s, L) = ((s << 2) & ((1 << L) - 1)) | (s >> (L - 2))

"""
    group_orbits(L, q)

Orbits of the configurations of `L` spins with `q` down spins under cell translations (and global
spin flip when `2q == L`), with each orbit's images and stabilizer; returns `(orbits, useflip)`.
The sparse 20-spin blocks need this physical-frame Sz basis: the core builds dense blocks only.
"""
function group_orbits(L, q)
    n = L ÷ 2; mask = (1 << L) - 1
    states = [s for s in 0:mask if count_ones(s) == q]
    seen = falses(1 << L); orbits = Any[]
    useflip = (2q == L)
    for rep in states
        seen[rep+1] && continue
        images = Dict{Int,Tuple{Int,Int}}(); stabilizers = Tuple{Int,Int}[]
        for f in 0:Int(useflip)
            s = f == 0 ? rep : xor(rep, mask)
            for t in 0:n-1
                if haskey(images, s)
                    t0, f0 = images[s]
                    push!(stabilizers, (t - t0, f - f0))
                else
                    images[s] = (t, f); seen[s+1] = true
                end
                s = rotate_cell(s, L)
            end
        end
        push!(orbits, (; rep, images, stabilizers))
    end
    @assert sum(length(o.images) for o in orbits) == binomial(L, q)
    return orbits, useflip
end

export sparse_block, low_spectrum, run_diagnostic, summarize_diagnostics

"Sparse Sz/cell-momentum/spin-flip block; numerical diagnostics only."
function sparse_block(delta, L, k, parity, orbits)
    n = L÷2
    odd = parity == -1 ? 1 : 0
    allowed = [o for o in orbits if
        all(mod(2k*t+n*odd*f, 2n) == 0 for (t,f) in o.stabilizers)]
    d = length(allowed)
    labels = zeros(Int, 1<<L)
    coeff = zeros(ComplexF64, 1<<L)
    for (i,o) in enumerate(allowed), (s,(t,f)) in o.images
        labels[s+1] = i
        coeff[s+1] = cispi(2k*t/n)*parity^f/sqrt(length(o.images))
    end
    rows, cols, values = Int[], Int[], ComplexF64[]
    for (j,o) in enumerate(allowed)
        s = o.rep
        scale = sqrt(length(o.images))
        diag = 0.0
        for site in 0:L-1
            next = mod(site+1, L)
            J = 1+(iseven(site) ? delta : -delta)
            opposite = ((s>>site)&1) != ((s>>next)&1)
            diag += J*(opposite ? -0.25 : 0.25)
            if opposite
                target = xor(s, (1<<site)|(1<<next))
                i = labels[target+1]
                i == 0 && continue
                push!(rows,i); push!(cols,j); push!(values,J/2*scale*conj(coeff[target+1]))
            end
        end
        push!(rows,j); push!(cols,j); push!(values,diag)
    end
    H = sparse(rows,cols,values,d,d)
    maximum(abs,H-H';init=0.0) < 1e-10 || error("non-Hermitian block")
    return Hermitian(H)
end

"Enumerate every sector, retaining numerical low Ritz values in large blocks."
function low_spectrum(delta,L; levels=40,dense_max=512)
    Random.seed!(20261014)
    BLAS.set_num_threads(1)
    records = Any[]
    E = Float64[]
    cutoff = Inf
    totaldim = 0
    for q in 0:L÷2
        orbits,flip = group_orbits(L,q)
        qdim = 0
        mult = flip ? 1 : 2
        conjugates = Dict{Tuple{Int,Int},Any}()
        for k in 0:L÷2-1, parity in (flip ? (1,-1) : (1,))
            start = time()
            if k > L÷4
                haskey(conjugates,(L÷2-k,parity)) || continue
                prior = conjugates[(L÷2-k,parity)]
                record = copy(prior)
                record["k"] = string(k)
                record["seconds"] = 0.0
                record["conjugate_of"] = string(L÷2-k)
                d = parse(Int,record["dimension"])
                qdim += d; totaldim += mult*d
                for _ in 1:mult
                    append!(E,record["ritz_values_float"])
                end
                push!(records,record)
                continue
            end
            H = sparse_block(delta,L,k,parity,orbits)
            d = size(H,1)
            d == 0 && continue
            qdim += d
            totaldim += mult*d
            complete = d <= dense_max
            residual = 0.0
            if complete
                vals = eigvals(Hermitian(Matrix(H)))
            else
                vals, vectors, info = eigsolve(H, randn(ComplexF64,d), levels, :SR;
                    ishermitian=true, krylovdim=max(128,3levels), maxiter=400, tol=1e-11)
                info.converged >= levels || error("Lanczos did not converge")
                vals = real.(vals[1:levels])
                residual = maximum(info.normres[1:levels])
                cutoff = min(cutoff, vals[end])
            end
            for _ in 1:mult
                append!(E,vals)
            end
            push!(records,Dict("q"=>string(q),"k"=>string(k),"parity"=>string(parity),
                "dimension"=>string(d),"multiplicity"=>string(mult),"complete"=>complete,
                "ritz_values_float"=>vals,"max_residual_float"=>residual,
                "seconds"=>time()-start))
            conjugates[(k,parity)] = records[end]
            println("diagnostic L=$L q=$q k=$k p=$parity d=$d sec=$(time()-start)")
            flush(stdout)
        end
        qdim == binomial(L,q) || error("charge dimension mismatch")
    end
    totaldim == 2^L || error("Hilbert dimension mismatch")
    # This heuristic tail assumes Ritz values include all lowest levels. It is NOT a proof.
    isfinite(cutoff) && filter!(e -> e <= cutoff, E)
    return PartialSpectrum(sort(E),2^L,isinf(cutoff)), records
end

"20-spin diagnostic with an explicitly ineligible odd seed n0=5."
function run_diagnostic(delta::Rational{BigInt})
    start = time()
    complete(cells) = float_ring_spectrum(rotated_model(delta), cells)
    small = complete(5)
    large,blocks = low_spectrum(Float64(delta),20)
    six, twelve = complete(3), complete(6)
    defect3 = twelve.E[1]-2six.E[1]
    defect5 = large.E[1]-2small.E[1]
    eight, sixteen = complete(4), complete(8)
    defect4 = sixteen.E[1]-2eight.E[1]
    xi34 = 1/log(defect3/defect4)
    xi45 = 1/log(defect4/defect5)
    rows = Any[]
    for bi in 1:128
        beta = bi/4
        zn = log_partition(small.E,beta)
        lo,hi = logZ_bounds(large,beta)
        l2,h2 = logZ_bounds(large,2beta)
        s = 1-exp(lo-2zn)
        t = 1-exp(l2-2hi)
        b = s < 1 && t < 1 ? certified_gap_bound(max(0.0,s),max(0.0,t),beta).bound : nothing
        push!(rows,Dict("beta_float"=>beta,"s_float"=>s,"t_float"=>t,"bound_float"=>b,
            "logZ_tail_width_float"=>hi-lo))
    end
    good = filter(r -> r["bound_float"] !== nothing,rows)
    root = normpath(joinpath(MODEL_DIR, "..", ".."))
    result = Dict("delta"=>string(delta),"spins"=>"20","n0_cells"=>"5",
        "theorem_eligible"=>false,"certified"=>false,
        "reason"=>"Odd seed; partial Float64 Ritz spectra can miss levels or multiplicities",
        "gap_20_float"=>large.E[2]-large.E[1],"sectors"=>blocks,
        "ground_energy_defects_float"=>[defect3,defect4,defect5],
        "decay_length_cells_3_4_float"=>xi34,"decay_length_cells_4_5_float"=>xi45,
        "decay_times_delta_two_thirds_float"=>xi45*Float64(delta)^(2/3),
        "scaling_estimate_n0_float"=>4*(0.25/Float64(delta))^(2/3),
        "scaling_caveat"=>"Energy-defect decay fit, not a certified correlation length or seed",
        "retained_levels"=>string(length(large.E)),"scan_float"=>rows,
        "best_odd_seed_float"=>isempty(good) ? nothing : argmax(r -> r["bound_float"],good),
        "seconds_total"=>time()-start,
        "metadata"=>run_metadata(Dict("delta"=>string(delta));
            source_sha256=source_hashes(root;dirs=(),files=(
                "models/dimerized_heisenberg/src/Diagnostics.jl",
                "models/dimerized_heisenberg/src/DimerizedHeisenberg.jl")),
            packages=(CertifiedGaps,KrylovKit),status="20_spin_diagnostic_only",certified=false))
    path = joinpath(@__DIR__, "..", "results",
        "diagnostic20_$(numerator(delta))_$(denominator(delta)).json")
    write_result(path,result)
    println("SAVED $path in $(time()-start) seconds")
end

"Scalar fields of a diagnostic run kept in the tracked summary (Ritz values and scan dropped)."
const SUMMARY_KEYS = ("delta", "spins", "n0_cells", "theorem_eligible", "certified", "reason",
    "gap_20_float", "ground_energy_defects_float", "decay_length_cells_3_4_float",
    "decay_length_cells_4_5_float", "decay_times_delta_two_thirds_float",
    "scaling_estimate_n0_float", "scaling_caveat", "retained_levels", "best_odd_seed_float",
    "seconds_total")

"""
    summarize_diagnostics(deltas)

Condense results/diagnostic20_<p>_<q>.json for every delta into results/diagnostic20_summary.json.
"""
function summarize_diagnostics(deltas)
    rows = Any[]
    for delta in deltas
        name = "$(numerator(delta))_$(denominator(delta))"
        full = parse_exact_json(read(joinpath(@__DIR__, "..", "results",
            "diagnostic20_$name.json"), String))
        row = Dict(k => full[k] for k in SUMMARY_KEYS)
        row["sectors"] = length(full["sectors"])
        row["complete_sectors"] = count(b -> b["complete"], full["sectors"])
        tofloat(x) = x isa AbstractString ? parse(Float64, x) : Float64(x)
        row["max_ritz_residual_float"] =
            maximum(tofloat(b["max_residual_float"]) for b in full["sectors"])
        row["run_metadata"] = full["metadata"]
        push!(rows, row)
    end
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    write_result(joinpath(@__DIR__, "..", "results", "diagnostic20_summary.json"), Dict(
        "scope" => "Numerical only: partial Float64 Ritz spectra of the 20-spin ring, an odd " *
            "n0 = 5 seed outside the even-seed theorem, energy-defect decay fits",
        "diagnostics" => rows,
        "metadata" => run_metadata(Dict("model" => "dimerized_heisenberg");
            source_sha256=source_hashes(root; dirs=(), files=(
                "models/dimerized_heisenberg/src/Diagnostics.jl",)),
            packages=(CertifiedGaps,), status="20_spin_diagnostic_summary", certified=false)))
end
end
