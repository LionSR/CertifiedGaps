"""Fresh, exact-input certificates and finite-grid rejection witnesses for the ladder."""
module LadderSuite

include("../certify.jl")
using CertifiedGaps.ExactArithmetic: round_up, round_down
using CertifiedGaps.RingModels

const Q = Rational{BigInt}
const GRID = big(10)^8
const BASE = normpath(joinpath(@__DIR__, ".."))
const BETAS = [big(k)//4 for k in 1:64]
const COUPLINGS = [4//1, 7//2, 3//1, 5//2, 2//1, 3//2,
    326133//131072, 163069//65536]
const FORMAT = "heisenberg-ladder verified suite v1"

q(x::AbstractString) = parse(Q, x)
q(x::Integer) = Q(x)
token(x) = string(Q(x))
endpoints(x) = [token(inf(x)), token(sup(x))]

"Exact rotated rung model and explicit real reflection-positive factors."
function rp_model(J)
    x = Q[0 1; 1 0]/2
    a = Q[0 1; -1 0]/2
    z = Q[1 0; 0 -1]/2
    id = Matrix{Q}(I, 2, 2)
    factors = vcat([kron(c, id) for c in (x,a,z)], [kron(id,c) for c in (x,a,z)])
    bond = -sum(kron(c,c) for c in factors)
    onsite = -J * (kron(x,x) + kron(a,a) + kron(z,z))
    model = RingModel(4, bond, onsite; name="checkerboard-rotated Heisenberg ladder")
    return model, rp_certificate(model; factors)
end

"Enclose both true seed defects using every verified eigenvalue and its multiplicity."
function defects(specs, n, beta)
    a = enclosed_logZ(specs[n], beta)
    b = enclosed_logZ(specs[2n], beta)
    c = enclosed_logZ(specs[2n], 2beta)
    return -expm1(b-2a), -expm1(c-2b), [a,b,c]
end

"Prove failure of attraction to zero from lower defects, or return no conclusion."
function rejection(s::Q, t::Q; steps=64)
    0 <= s <= 1 && 0 <= t <= 1 || error("invalid lower defects")
    grid = big(10)^24
    for j in 0:steps
        z = 1 - (1-t)*(1-s)^2
        # z >= 1/2 implies both next components >= 1/2, hence no attraction to zero.
        z >= 1//2 && return Dict("step"=>token(j), "z_lower"=>token(z))
        sp = round_down(z^2/(2*(1-z)^2), grid)
        zt = 1 - (1-sp)*(1-t)^2
        zt >= 1//2 && return Dict("step"=>token(j), "thermal_z_lower"=>token(zt))
        s, t = sp, round_down(zt^2/(2*(1-zt)^2), grid)
    end
    return nothing
end

"Classify one point rigorously; a failed upper-bound attempt is never called rejection."
function classify(specs, n, beta)
    s, t, logs = defects(specs, n, beta)
    su = max(1//GRID, round_up(Q(sup(s)), GRID))
    tu = max(1//GRID, round_up(Q(sup(t)), GRID))
    sl = max(big(0)//1, round_down(Q(inf(s)), GRID))
    tl = max(big(0)//1, round_down(Q(inf(t)), GRID))
    answer = su < 1 && tu < 1 ? certified_gap_bound(su,tu,beta).bound : nothing
    witness = rejection(sl,tl)
    answer !== nothing && witness !== nothing && error("inconsistent classification")
    status = answer !== nothing ? "certified" :
        witness !== nothing ? "rejected_by_recurrence" : "unresolved"
    record = Dict{String,Any}("n0"=>token(n), "beta0"=>token(beta), "status"=>status,
        "eps_s_upper"=>token(su), "eps_t_upper"=>token(tu),
        "eps_s_lower"=>token(sl), "eps_t_lower"=>token(tl))
    if answer !== nothing
        record["gap_lower_bound"] = token(round_down(Q(answer), big(10)^9))
        record["interval_bound_exact"] = token(answer)
        record["logZ"] = endpoints.(logs)
    end
    witness !== nothing && (record["rejection"] = witness)
    return record
end

"Outward seed continuation on an entire coupling interval using centered rung norms."
function continuation(seed, delta)
    n, beta = Int(q(seed["n0"])), q(seed["beta0"])
    s, t = q(seed["eps_s_upper"]), q(seed["eps_t_upper"])
    su = round_up(Q(sup(1-(1-to_interval(s))*exp(-to_interval(2beta*n*delta)))), GRID)
    tu = round_up(Q(sup(1-(1-to_interval(t))*exp(-to_interval(4beta*n*delta)))), GRID)
    b = certified_gap_bound(su,tu,beta).bound
    b === nothing && return nothing
    return Dict("delta"=>token(delta), "eps_s_upper"=>token(su), "eps_t_upper"=>token(tu),
        "gap_lower_bound"=>token(round_down(Q(b),big(10)^9)))
end

"Ground energy bounds from the complete spectrum, with no ordering assumption on radii."
ground_bounds(s) = (minimum(Q(e)-Q(r) for (e,r) in zip(s.energies,s.radii)),
    minimum(Q(e)+Q(r) for (e,r) in zip(s.energies,s.radii)))

"Cache rigorously shifted trace endpoints for continuous-temperature exclusion."
function temperature_data(specs,n)
    a,b = ground_bounds(specs[n]),ground_bounds(specs[2n])
    cache = Dict{Tuple{Int,Q},Interval{Float64}}()
    function logbar(L,beta)
        get!(cache,(L,beta)) do
            shift = L == n ? a[1] : b[1]
            enclosed_logZ(specs[L],beta) + to_interval(beta*shift)
        end
    end
    return (; a,b,logbar)
end

"Lower defects valid on [lo,hi]; spatial monotonicity is not assumed."
function window_lower(data,n,lo,hi)
    # H_L - shift_L is PSD, so each shifted partition sum decreases with beta.
    beta = interval(Float64(lo),Float64(hi))
    @assert Q(Float64(lo)) == lo && Q(Float64(hi)) == hi
    logps = data.logbar(2n,lo) - 2data.logbar(n,hi) - beta*to_interval(data.b[1]-2data.a[1])
    # Thermal purity increases with beta: d log P_tau / d beta = 2(E_beta-E_2beta)>=0.
    logpt = data.logbar(2n,2hi) - 2data.logbar(2n,hi)
    sl = max(big(0)//1,round_down(Q(inf(-expm1(logps))),GRID))
    tl = max(big(0)//1,round_down(Q(inf(-expm1(logpt))),GRID))
    return sl,tl
end

"Reject all beta >= B with a positive ground-energy mismatch and a thermal trace bound."
function tail_lower(data,n,B)
    delta = data.b[1]-2data.a[2]
    delta > 0 || return nothing
    logps = data.logbar(2n,B)-to_interval(B*delta)
    sl = max(big(0)//1,round_down(Q(inf(-expm1(logps))),GRID))
    rejection(sl,big(0)//1) === nothing && return nothing
    return Dict("beta_from"=>token(B), "ground_mismatch_lower"=>token(delta),
        "eps_s_lower"=>token(sl), "eps_t_lower"=>"0//1")
end

"Prove that neither eligible seed works at any positive beta, or fail explicitly."
function exclude_temperatures(specs,n)
    data = temperature_data(specs,n)
    B = big(1)//1
    tail = tail_lower(data,n,B)
    while tail === nothing && B < 1024
        B *= 2
        tail = tail_lower(data,n,B)
    end
    tail === nothing && error("could not exclude the infinite-temperature-parameter tail")
    windows = Any[]
    function split(lo,hi,depth)
        sl,tl = window_lower(data,n,lo,hi)
        witness = rejection(sl,tl)
        if witness !== nothing
            push!(windows,Dict("beta_from"=>token(lo),"beta_to"=>token(hi),
                "eps_s_lower"=>token(sl),"eps_t_lower"=>token(tl),"rejection"=>witness))
        else
            depth < 20 || error("temperature window unresolved: $lo to $hi")
            split(lo,(lo+hi)/2,depth+1)
            split((lo+hi)/2,hi,depth+1)
        end
    end
    split(big(0)//1,B,0)
    return Dict("n0"=>token(n),"windows"=>windows,"tail"=>tail)
end

"Check interval coverage and every lower-defect witness including the infinite-beta tail."
function check_exclusion(specs,record)
    n = Int(q(record["n0"]))
    n in (2,4) || error("invalid excluded seed length")
    data = temperature_data(specs,n)
    previous = big(0)//1
    for row in record["windows"]
        lo,hi = q(row["beta_from"]),q(row["beta_to"])
        lo == previous && hi > lo || error("temperature coverage gap")
        sl,tl = window_lower(data,n,lo,hi)
        s,t = q(row["eps_s_lower"]),q(row["eps_t_lower"])
        s <= sl && t <= tl || error("temperature lower defect fails")
        rejection(s,t) == row["rejection"] || error("temperature rejection fails")
        previous = hi
    end
    tail = record["tail"]
    previous == q(tail["beta_from"]) || error("tail coverage gap")
    fresh = tail_lower(data,n,previous)
    fresh !== nothing || error("tail rejection fails")
    q(tail["ground_mismatch_lower"]) <= q(fresh["ground_mismatch_lower"]) ||
        error("ground mismatch fails")
    q(tail["eps_s_lower"]) <= q(fresh["eps_s_lower"]) || error("tail defect fails")
    q(tail["eps_t_lower"]) == 0 || error("invalid tail thermal defect")
    rejection(q(tail["eps_s_lower"]),big(0)//1) !== nothing || error("tail witness fails")
    return true
end

"Source and environment fingerprints for the current model project and shared core."
function metadata(config)
    hashes = source_hashes(ROOT; dirs=("julia/src", "models/heisenberg_ladder"),
        files=("julia/Project.toml", "models/heisenberg_ladder/Project.toml",
            "models/heisenberg_ladder/Manifest.toml", "outputs/thermal-mpo-gap-scope.tex"))
    return run_metadata(config; source_sha256=hashes,
        packages=(CertifiedGaps,IntervalArithmetic,JSON3), status="certified", certified=true)
end

"Compute the exact orbit block dimensions, independently of any saved spectrum file."
function sector_dimensions(L)
    rows = Dict{NTuple{4,Int},Tuple{Int,Int}}()
    for m in 0:L
        o = orbit_data(L,m)
        for k in 0:fld(L,2), p in 0:1, r in 0:Int(o.flip)
            selected, _ = block_basis(o,k,p,r)
            isempty(selected) && continue
            mult = (m==L ? 1 : 2)*(mod(2k,L)==0 ? 1 : 2)
            rows[(m,k,(-1)^p,o.flip ? (-1)^r : 0)] = (length(selected),mult)
        end
    end
    @assert sum(d*m for (d,m) in values(rows)) == 4^L
    return rows
end

"Recompute one coupling; the output includes every beta-grid decision, including failures."
function certify_coupling(J)
    start = time()
    @assert exact_rp(J)
    model, rp = rp_model(J)
    specs = Dict(L=>full_spectrum(L,J; verified=true) for L in (2,4,8))
    rows = [classify(specs,n,beta) for n in (2,4) for beta in BETAS]
    good = filter(r->r["status"]=="certified", rows)
    best = isempty(good) ? nothing : argmax(r->q(r["gap_lower_bound"]), good)
    cont = best === nothing ? nothing : continuation(best, big(1)//10000)
    excluded = J in (2//1,3//2) ? [exclude_temperatures(specs,n) for n in (2,4)] : nothing
    summaries = Dict(string(L)=>Dict("dimension"=>token(s.dimension),
        "max_radius"=>token(round_up(Q(maximum(s.radii)),GRID)), "seconds"=>s.seconds,
        "block_count"=>token(length(s.rows))) for (L,s) in specs)
    return Dict("J"=>token(J), "model"=>encode_model(model), "rp"=>encode_rp(rp),
        "grid"=>rows, "best"=>best, "continuation"=>cont, "spectra"=>summaries,
        "all_beta_exclusions"=>excluded,
        "seconds"=>time()-start)
end

"Produce a fresh suite without trusting any numerical scan or old certificate as input."
function produce(path=joinpath(BASE,"results/certified_suite.json"))
    start = time()
    records = Any[]
    for J in COUPLINGS
        result = certify_coupling(J)
        push!(records,result)
        println("COUPLING ",J," best=",result["best"] === nothing ? "none" :
            result["best"]["gap_lower_bound"], " seconds=",result["seconds"])
        flush(stdout)
    end
    cert = Dict("format"=>FORMAT, "parameters"=>records,
        "metadata"=>metadata(Dict("n0"=>["2//1","4//1"],
            "beta_grid"=>"1//4:1//4:16//1", "source_of_spectra"=>"fresh complete verification")),
        "validity"=>"even L >= n0 rungs; N=2L spin sites, N divisible by 4",
        "seconds"=>time()-start)
    write_result(path,cert)
    return cert
end

"Check an entire certificate from fresh spectra, exact RP, and the certified recurrence."
function check(path=joinpath(BASE,"results/certified_suite.json"))
    start = time()
    cert = parse_exact_json(read(path,String))
    cert["format"] == FORMAT || error("unknown certificate format")
    cert["validity"] == "even L >= n0 rungs; N=2L spin sites, N divisible by 4" ||
        error("incorrect validity claim")
    for (file,hash) in cert["metadata"]["source_sha256"]
        sha256_file(joinpath(ROOT,file)) == hash || error("source mismatch: $file")
    end
    q.(getindex.(cert["parameters"],"J")) == Q.(COUPLINGS) || error("couplings mismatch")
    for result in cert["parameters"]
        J = Rational{Int}(q(result["J"]))
        @assert exact_rp(J)
        model,rp = rp_model(J)
        stored = decode_model(result["model"])
        stored.d == model.d && stored.bond == model.bond && stored.onsite == model.onsite ||
            error("model mismatch")
        verify_rp(stored,decode_rp(result["rp"])) || error("RP identity fails")
        specs = Dict(L=>full_spectrum(L,J; verified=true) for L in (2,4,8))
        expected = [(n,beta) for n in (2,4) for beta in BETAS]
        length(result["grid"]) == length(expected) || error("incomplete grid")
        for (row,(n,beta)) in zip(result["grid"],expected)
            (Int(q(row["n0"])),q(row["beta0"])) == (n,beta) || error("grid mismatch")
            s,t,logs = defects(specs,n,beta)
            sl,tl = q(row["eps_s_lower"]),q(row["eps_t_lower"])
            su,tu = q(row["eps_s_upper"]),q(row["eps_t_upper"])
            sl <= max(0,Q(inf(s))) && tl <= max(0,Q(inf(t))) || error("lower defects fail")
            Q(sup(s)) <= su && Q(sup(t)) <= tu || error("upper defects fail")
            if row["status"] == "certified"
                b = certified_gap_bound(su,tu,beta).bound
                b !== nothing && 0 < q(row["gap_lower_bound"]) <= Q(b) ||
                    error("gap bound fails")
                q(row["interval_bound_exact"]) <= Q(b) || error("interval bound mismatch")
                for (saved,fresh) in zip(row["logZ"],logs)
                    # Independent eigensolves need overlap; the seeds are checked above.
                    q(saved[1]) <= Q(sup(fresh)) && Q(inf(fresh)) <= q(saved[2]) ||
                        error("log Z enclosure mismatch")
                end
            elseif row["status"] == "rejected_by_recurrence"
                witness = rejection(sl,tl)
                witness !== nothing || error("rejection not proved")
                witness == row["rejection"] || error("rejection witness mismatch")
            else
                row["status"] == "unresolved" || error("unknown status")
            end
        end
        good = filter(r->r["status"]=="certified", result["grid"])
        best = isempty(good) ? nothing : argmax(r->q(r["gap_lower_bound"]),good)
        best == result["best"] || error("best seed mismatch")
        if result["continuation"] !== nothing
            fresh = continuation(best,q(result["continuation"]["delta"]))
            fresh == result["continuation"] || error("continuation mismatch")
        end
        if J in (2//1,3//2)
            exclusions = result["all_beta_exclusions"]
            q.(getindex.(exclusions,"n0")) == [2,4] || error("exclusion lengths mismatch")
            foreach(r->check_exclusion(specs,r),exclusions)
        else
            result["all_beta_exclusions"] === nothing || error("unexpected exclusion")
        end
        for (L,s) in specs
            summary = result["spectra"][string(L)]
            Int(q(summary["dimension"])) == s.dimension == 4^L || error("dimension mismatch")
            Int(q(summary["block_count"])) == length(s.rows) || error("block count mismatch")
            Q(maximum(s.radii)) <= q(summary["max_radius"]) || error("spectral radius mismatch")
        end
    end
    println("PASS: fresh RP, complete spectra, seed bounds, grid decisions and continuation; ",
        round(time()-start; digits=2)," seconds")
    return true
end

end
