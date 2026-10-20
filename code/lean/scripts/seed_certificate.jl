# Rational data for `CertifiedGaps.TableCertificate` (lean/CertifiedGaps/Seed/Certificate.lean): a
# table of twisted ring traces split into the sectors of a `SectorScheme` (plain, D₂, or I/O/N).
#
# Given a trace table (centers c_g(n) for every twist g of the scheme), the trace error ε that
# the filters accept, the seed shape (moment order q, seed length n₀ = q a, dominance threshold
# ℓ) and the parameters (β, κ), this searches one polynomial filter X^(2a') q(X)² per capped
# sector in floating point, then computes in exact rational arithmetic the filter bounds, the caps (the
# smallest `digits`-digit decimals passing the ray test), the masses and both seed defects, and
# prints the Lean certificate. Lean re-checks every inequality by `decide`; this script only
# searches.
#
# Usage: julia lean/scripts/seed_certificate.jl config.toml
#
# The config (TOML; every rational is a string literal parsed exactly, never evaluated):
#   name = "seed36"             # Lean name of the certificate
#   scheme = "d2"               # "plain", "d2" or "ion" (SectorScheme in Seed/Scheme.lean)
#   table_json = "models/spin1/results/source1.json"   # rows with "n", "twist", "center", or
#   [center.0] / [center.1] ... # inline tables "4" = "124885379/1000000", ... per twist
#   lengths = [4, 5, 6, 7, 8, 9, 10, 11, 12]
#   error = "1/20000"           # trace error accepted by the filters (> the certified error)
#   q = 12                      # moment order of the masses (even, a table length)
#   a = 3                       # seed length n₀ = q a
#   ell = "999/1000"            # dominance threshold
#   beta = "49/4"
#   kappa = "700741/500000"
#   spatial_defect = "47916/1000000"   # optional: default = achieved value rounded up
#   thermal_defect = "181521/1000000"  # optional
#   capped = [true, true]       # optional: which sectors get a filter (sector 0 must)
#   odd = [true, false]         # optional: which sector filters may use odd lengths
#   digits = 10                 # significant digits of the caps
# Paths are relative to the repository root, the working directory of the command.

include(joinpath(@__DIR__, "chain_certificate.jl"))  # parse_rational, round_up, round_down
using LinearAlgebra, TOML

const Q = Rational{BigInt}

"Read the table: `Dict(twist => Dict(n => center))` from `table_json` or inline `[center.<g>]`."
function read_table(cfg)
    tables = Dict{Int,Dict{Int,Q}}()
    if haskey(cfg, "table_json")
        # Minimal reader for rows {"n": .., "twist": .., "center": ".."}: no JSON dependency.
        text = read(cfg["table_json"], String)
        for m in eachmatch(r"\{[^{}]*\}", text)
            row = m.match
            n = match(r"\"n\"\s*:\s*\"?(\d+)", row)
            tw = match(r"\"twist\"\s*:\s*\"?(\d+)", row)
            ce = match(r"\"center\"\s*:\s*\"([^\"]+)\"", row)
            (n === nothing || tw === nothing || ce === nothing) && continue
            get!(tables, parse(Int, tw[1]), Dict{Int,Q}())[parse(Int, n[1])] =
                parse_rational(ce[1])
        end
    else
        for (g, tab) in cfg["center"], (k, v) in tab
            get!(tables, parse(Int, g), Dict{Int,Q}())[parse(Int, k)] = parse_rational(v)
        end
    end
    return tables
end

"Coefficients (constant term first) of `X^shift q(X)²`, as Lean's `filterPoly (shift/2) q`."
function convsquare(q, shift)
    out = zeros(eltype(q), 2length(q) - 1 + shift)
    for i in eachindex(q), j in eachindex(q)
        out[i+j-1+shift] += q[i] * q[j]
    end
    return out
end

"Taylor coefficients of `p(±(v + y))` in `y`, constant term first."
taylor(P, v, sign) = [sum(binomial(i, j) * P[i+1] * sign^i * v^(i - j) for i in j:length(P)-1)
    for j in 0:length(P)-1]

"Lean's `CapCheck`: on both rays the Taylor coefficients are ≥ 0 and the constant term > B."
cap_ok(P, B, v) = all(s -> (t = taylor(P, v, s); t[1] > B && all(>=(0), t[2:end])), (-1, 1))

"Lean's `pairingBound` for the symmetric intervals `m_k ± e`."
pairing_bound(P, m, e) = sum((c > 0 ? c * (m[k-1] + e) : c * (m[k-1] - e)) for (k, c) in
    enumerate(P) if c != 0; init = zero(Q))

"Floating-point search for `X^s q²` (s ≥ 4 even): least-squares `q` on the moment Hankel matrix."
function search_filter(mom::Dict{Int,Float64}, eps, nmax; allowodd = true)
    best = nothing
    for shift in 4:2:nmax-2, stride in (allowodd ? (1, 2) : (2,))
        powers = collect(0:stride:div(nmax - shift, 2))
        length(powers) >= 2 || continue
        all(shift + i + j in keys(mom) for i in powers for j in powers) || continue
        M = [mom[shift+i+j] for i in powers, j in powers]
        for ridge in (0.0, eps, 3eps, 10eps, 100eps), at in (1.005, 1.0, 0.95, 0.9, 0.875, 0.85, 0.8)
            A = Symmetric(M + ridge * I)
            minimum(eigvals(A)) > 1e-13 || continue
            qq = A \ (at .^ powers)
            qq ./= maximum(abs, qq)
            qv = zeros(maximum(powers) + 1)
            qv[powers.+1] = qq
            P = convsquare(qv, shift)
            B = sum(c * mom[i-1] for (i, c) in enumerate(P) if c != 0) + eps * sum(abs, P)
            B > 0 || continue
            lo, hi = 0.05, 2.0
            cap_ok(P, B, hi) || continue
            for _ in 1:50
                mid = (lo + hi) / 2
                cap_ok(P, B, mid) ? (hi = mid) : (lo = mid)
            end
            (best === nothing || hi < best.cap) && (best = (cap = hi, shift, q = qv))
        end
    end
    best === nothing && error("no filter found")
    return best
end

"Exact filter: rationalized `q`, its bound, and the smallest `digits`-digit cap."
function exact_filter(fq, shift, m::Dict{Int,Q}, e::Q, digits)
    q = [rationalize(BigInt, x; tol = 1e-9) for x in fq.q]
    P = convsquare(q, shift)
    B = pairing_bound(P, m, e)
    lo, hi = big(1) // 20, big(2) // 1
    cap_ok(P, B, hi) || error("rationalized filter fails at v = 2")
    for _ in 1:200
        mid = (lo + hi) / 2
        cap_ok(P, B, mid) ? (hi = mid) : (lo = mid)
    end
    cap = round_up(hi, digits)
    cap_ok(P, B, cap) || error("rounded cap fails the ray test")
    return (q = q, a = div(shift, 2), B = B, cap = cap)
end

"Lean literal of a rational: `m` or `m / n` (parenthesized when negative inside a list)."
leanq(x::Q) = denominator(x) == 1 ? string(numerator(x)) : string(numerator(x), " / ", denominator(x))
leanlist(xs) = "[" * join(leanq.(xs), ", ") * "]"

function lean_table(name, c)
    lines = ["def $name : ℕ → ℚ"]
    for n in sort(collect(keys(c)))
        push!(lines, "  | $n => $(leanq(c[n]))")
    end
    push!(lines, "  | _ => 0")
    return join(lines, "\n")
end

"The sector schemes of `lean/CertifiedGaps/Seed/Scheme.lean`: `Z_g = ∑_s coeff[g, s] M_s`."
const SCHEMES = Dict(
    "plain" => (coeff = [1;;], inv = Q[1;;]),
    "d2" => (coeff = [1 3; 1 -1], inv = Q[1//4 3//4; 1//4 -1//4]),
    "ion" => (coeff = [1 2 3; 1 2 -1; 1 -1 0],
        inv = Q[1//12 3//12 8//12; 1//12 3//12 -4//12; 1//4 -1//4 0]))

"Lean's `SectorSeed.sectorBound`: the `q r`th moment bound of a nontrivial sector."
function sector_bound(cap, mass, q, r)
    cap === nothing && return mass^r
    h = cap^q
    return mass <= h ? mass^r : h^r + (mass - h)^r
end

function main(path)
    cfg = TOML.parsefile(path)
    name = cfg["name"]
    scheme = SCHEMES[get(cfg, "scheme", "d2")]
    nt, ns = size(scheme.coeff)
    tables = read_table(cfg)
    lengths = Int.(cfg["lengths"])
    ε = parse_rational(cfg["error"])
    q, a = cfg["q"], cfg["a"]
    ℓ = parse_rational(cfg["ell"])
    digits = get(cfg, "digits", 10)
    capped = Bool.(get(cfg, "capped", fill(true, ns)))
    odd = Bool.(get(cfg, "odd", [s == 1 for s in 1:ns]))
    capped[1] || error("the dominant sector needs a filter")
    nmax = maximum(lengths)
    # Sector moments M_s = ∑_g inv[s, g] Z_g and their errors ε ∑_g |inv[s, g]|, exactly as the
    # interval bounds `loS`, `hiS` of the symmetric table in Lean.
    M = [Dict(n => sum(scheme.inv[s, g] * tables[g-1][n] for g in 1:nt) for n in lengths)
        for s in 1:ns]
    e = [ε * sum(abs, scheme.inv[s, :]) for s in 1:ns]
    filters = Any[nothing for _ in 1:ns]
    for s in 1:ns
        capped[s] || continue
        mom = Dict(n => Float64(v) for (n, v) in M[s] if odd[s] || iseven(n))
        f = search_filter(mom, Float64(e[s]), nmax; allowodd = odd[s])
        filters[s] = exact_filter(f, f.shift, M[s], e[s], digits)
    end
    mult = scheme.coeff[1, :]
    caps = [f === nothing ? nothing : f.cap for f in filters]
    mass = [M[s][q] + e[s] for s in 1:ns]
    hJ = ℓ^q
    @assert hJ <= mass[1] <= 2hJ "the dominant mass must lie in [ℓ^q, 2ℓ^q]"
    for s in 2:ns
        caps[s] === nothing || mass[s] <= 2caps[s]^q || error("sector $s: mass > 2 cap^q")
    end
    R(r) = sum((mult[s] * sector_bound(caps[s], mass[s], q, r) for s in 2:ns); init = zero(Q)) +
        (mass[1] - hJ)^r
    dom = hJ^(2a) + R(2a)
    dom < 1 || error("no dominant eigenvalue: capped masses give Z_{2n₀} ≤ $(Float64(dom))")
    p0 = 1 - 1 / (1 + R(a) / hJ^a)^2
    thermal = (caps[1]^q)^(2a) + R(2a)
    q0 = 1 - 1 / thermal^2
    p = haskey(cfg, "spatial_defect") ? parse_rational(cfg["spatial_defect"]) : round_up(p0, 6)
    t = haskey(cfg, "thermal_defect") ? parse_rational(cfg["thermal_defect"]) : round_up(q0, 6)
    (1 - p) * (1 + R(a) / hJ^a)^2 <= 1 || error("spatial defect $(Float64(p)) < $(Float64(p0))")
    (1 - t) * thermal^2 <= 1 || error("thermal defect $(Float64(t)) < $(Float64(q0))")
    println("-- achieved defects: p₀ = ", Float64(p0), ", q₀ = ", Float64(q0),
        "; dominance bound ", Float64(dom))
    for g in 1:nt
        println(lean_table(name * "Center$(g-1)", Dict(n => tables[g-1][n] for n in lengths)), "\n")
    end
    opt(x) = x === nothing ? "none" : "some ($(leanq(x)))"
    optf(f) = f === nothing ? "none" : "some ($(f.a), $(leanlist(f.q)))"
    vec(xs) = "![" * join(xs, ", ") * "]"
    out = [
        "def $name : TableCertificate where",
        "  scheme := SectorScheme.$(get(cfg, "scheme", "d2"))",
        "  table := TraceBounds.centered $(vec(["$(name)Center$(g-1)" for g in 1:nt]))",
        "    $(lengths) ($(leanq(ε)))",
        "  filter := ![" * join(optf.(filters), ",\n    ") * "]",
        "  seed :=",
        "    { q := $q, a := $a, mult := $(vec(mult))",
        "      cap := $(vec(opt.(caps)))",
        "      mass := $(vec(leanq.(mass)))",
        "      ell := $(leanq(ℓ))",
        "      spatialDefect := $(leanq(p)) }",
        "  beta := $(leanq(parse_rational(cfg["beta"])))",
        "  kappa := $(leanq(parse_rational(cfg["kappa"])))",
        "  thermalDefect := $(leanq(t))",
        "",
        "-- Each condition is checked by its own `decide`: one `decide` on the whole conjunction",
        "-- exceeds the kernel's recursion depth.",
        "theorem $(name)_valid : $(name).Valid :=",
        "  ⟨by decide +kernel, by decide +kernel, by decide +kernel, by decide +kernel,",
        "    by decide +kernel, by decide, by decide, by decide +kernel, by decide +kernel,",
        "    by decide +kernel⟩"]
    println(join(out, "\n"))
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main(ARGS[1])
end
