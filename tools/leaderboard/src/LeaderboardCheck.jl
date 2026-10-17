"""
    LeaderboardCheck

Validator for the leaderboard data `home_page/_data/records.yml` (schema in its header) and the
news list `home_page/_data/news.yml`: field formats, referenced paths, every displayed decimal
against the exact bound (logarithms enclosed by the outward-rounded rational `log_enclosure` of
`CertifiedGaps`), the tier of each headline bound against its Lean theorem, agreement with the
generated `lean_records.json`, and agreement of each row with the certificate it links.
"""
module LeaderboardCheck

using Dates
using JSON3
using YAML
using CertifiedGaps.ExactArithmetic: log_enclosure

include("certificates.jl")

export Report, check_records, check_news, parse_exact, parse_decimal, bound_enclosure, main

const TIERS = ("lean-unconditional", "lean-conditional", "certified", "target")
const LEAN_KIND = Dict("lean-unconditional" => "unconditional", "lean-conditional" => "conditional")
const ID = r"^[a-z0-9]+(-[a-zA-Z0-9]+)*$"
const LEAN_NAME = r"^Leaderboard\.Records(\.[A-Za-z_][A-Za-z0-9_']*)+$"
const HANDLE = r"^[A-Za-z0-9](?:[A-Za-z0-9-]{0,38})$"   # a GitHub login
const URL = r"^https?://[^\s\"<>]+$"

"Errors (which fail the check), warnings and notes collected while validating."
struct Report
    errors::Vector{String}
    warnings::Vector{String}
    notes::Vector{String}
end
Report() = Report(String[], String[], String[])
err!(r::Report, msg) = (push!(r.errors, msg); r)
passed(r::Report) = isempty(r.errors)

"""
    parse_exact(x) -> Rational{BigInt}

Exact rational from an integer or a string: `"13480"`, `"-3/4"`, `"97//4"`, a decimal `"413.4"`
or scientific notation `"4134e-1"` (as Lean prints `ℚ` literals). Floats are refused.
"""
parse_exact(x::Integer) = BigInt(x) // big(1)
function parse_exact(x::AbstractString)
    s = replace(strip(x), " " => "", "(" => "", ")" => "")
    m = match(r"^([+-]?\d+)//?(\d+)$", s)
    if m !== nothing
        den = parse(BigInt, m[2])
        den == 0 && throw(ArgumentError("zero denominator in \"$x\""))
        return parse(BigInt, m[1]) // den
    end
    m = match(r"^([+-]?)(\d+)(?:\.(\d+))?(?:[eE]([+-]?\d+))?$", s)
    m === nothing && throw(ArgumentError("not an exact number: \"$x\""))
    digits = something(m[3], "")
    q = parse(BigInt, m[2] * digits) // big(10)^length(digits)
    e = m[4] === nothing ? 0 : parse(Int, m[4])
    q = e >= 0 ? q * big(10)^e : q / big(10)^(-e)
    return m[1] == "-" ? -q : q
end
parse_exact(x) = throw(ArgumentError("not an exact number: $(repr(x)); write it as a string"))

"A displayed decimal such as `\"0.3922\"` as `(exact value, number of decimals)`."
function parse_decimal(s)
    m = s isa AbstractString ? match(r"^(\d+)\.(\d+)$", s) : nothing
    m === nothing && throw(ArgumentError("decimal $(repr(s)) must be a quoted string like \"0.1234\""))
    return parse(BigInt, m[1] * m[2]) // big(10)^length(m[2]), length(m[2])
end

"""
    bound_enclosure(bound) -> (lo, hi)

Rational enclosure of an exact bound: log(r)/β for `{r, beta}` (outward-rounded), or the point
`value`.
"""
function bound_enclosure(b::AbstractDict)
    if haskey(b, "value")
        (haskey(b, "r") || haskey(b, "beta")) && throw(ArgumentError("give {r, beta} or {value}, not both"))
        v = parse_exact(b["value"])
        v > 0 || throw(ArgumentError("value must be positive"))
        return v, v
    end
    (haskey(b, "r") && haskey(b, "beta")) || throw(ArgumentError("need {r, beta} or {value}"))
    r, β = parse_exact(b["r"]), parse_exact(b["beta"])
    (r > 1 && β > 0) || throw(ArgumentError("need r > 1 and beta > 0"))
    lo, hi = log_enclosure(r)
    return lo / β, hi / β
end

"""
    check_decimal!(report, where, decimal, lo, hi)

The decimal must be rounded down from the exact value in [lo, hi] (decimal ≤ lo) and be its floor
at the displayed precision (decimal + 10⁻ᵏ > hi).
"""
function check_decimal!(report, where, decimal, lo, hi)
    q, k = try
        parse_decimal(decimal)
    catch e
        return err!(report, "$where: $(sprint(showerror, e))")
    end
    q <= lo || err!(report, "$where: decimal $decimal exceeds the exact bound ≈ $(Float64(lo))")
    q + big(1) // big(10)^k > hi ||
        err!(report, "$where: decimal $decimal is not rounded at its precision (exact ≈ $(Float64(hi)))")
    return report
end

getk(d, k) = d isa AbstractDict ? get(d, k, nothing) : nothing
isstr(x) = x isa AbstractString && !isempty(strip(x))
need_str!(report, where, d, keys...) =
    foreach(k -> isstr(getk(d, k)) || err!(report, "$where: missing $k"), keys)

function check_path!(report, where, root, path; dir = false)
    path isa AbstractString || return err!(report, "$where: path must be a string")
    (isabspath(path) || occursin("..", path)) &&
        return err!(report, "$where: path \"$path\" must be relative to the repository root")
    (dir ? isdir : isfile)(joinpath(root, path)) || err!(report, "$where: \"$path\" does not exist")
    return report
end

check_link!(report, where, root, link) =
    occursin(URL, string(link)) ? report : check_path!(report, where, root, link)

function check_date!(report, where, d)
    d isa Date && return report
    d isa AbstractString && tryparse(Date, d, dateformat"yyyy-mm-dd") !== nothing && return report
    return err!(report, "$where: date $(repr(d)) is not yyyy-mm-dd")
end

function check_exact!(report, where, x; positive = true)
    try
        q = parse_exact(x)
        positive && q <= 0 && err!(report, "$where: must be positive")
    catch e
        err!(report, "$where: $(sprint(showerror, e))")
    end
    return report
end

"Lean records keyed by name, from the JSON array written by `lake exe leaderboard_audit`."
function load_lean_records(path)
    data = JSON3.read(read(path, String))
    data isa JSON3.Array || throw(ArgumentError("$path: expected a JSON array of records"))
    out = Dict{String,Any}()
    for rec in data
        name = String(rec[:name])
        haskey(out, name) && throw(ArgumentError("$path: duplicate Lean record $name"))
        out[name] = rec
    end
    return out
end

"Compare one Lean entry of records.yml with its record in lean_records.json; returns the kind."
function check_lean_json!(report, where, lean, L0, json)
    rec = get(json, lean["name"], nothing)
    rec === nothing && return err!(report, "$where: Lean record $(lean["name"]) is not in lean_records.json")
    get(rec, :L0, nothing) == L0 || err!(report, "$where: L0 = $L0 but $(get(rec, :L0, nothing)) in Lean")
    for key in ("r", "beta")
        try
            parse_exact(lean[key]) == parse_exact(string(rec[Symbol(key)])) ||
                err!(report, "$where: $key = $(lean[key]) but $(rec[Symbol(key)]) in Lean")
        catch e
            err!(report, "$where: cannot compare $key: $(sprint(showerror, e))")
        end
    end
    kind = get(rec, :kind, nothing)
    kind in ("unconditional", "conditional") || err!(report, "$where: unknown Lean kind $(repr(kind))")
    return kind
end

"""
    check_lean_seeds!(report, where, rec, lean, leanjson, rowcert, root)

A conditional Lean record assumes `P_s(n0, β) ≥ 1 - s` and `P_τ(2n0, β) ≥ 1 - t`. The certificate
must prove them: the same `n0` and `β`, and defects `εₛ ≤ s` and `ε_τ ≤ t`. The certificate is
the row's, or `lean.certificate` when the Lean theorem uses another seed. A certificate entry
without defects falls back to the one-step entry of the seed (`window_defects`). A two-seed
record (`seed2`) must find its second seed among the one-step entries of the same certificate.
Skipped while the Lean record does not state `n0`, `s` and `t`.
"""
function check_lean_seeds!(report, where, rec, lean, leanjson, rowcert, root)
    all(k -> haskey(leanjson, k), (:n0, :s, :t)) || return report
    c = if haskey(lean, "certificate")
        check_certificate!(report, "$where: lean", rec, nothing, root; path = lean["certificate"],
            seed = Dict("n0" => leanjson[:n0], "beta0" => lean["beta"]))
    else
        rowcert
    end
    c === nothing && return err!(report, "$where: no certificate entry for the Lean seed")
    (cert_exact(c.n0) == parse_exact(leanjson[:n0]) && cert_exact(c.beta) == parse_exact(lean["beta"])) ||
        return err!(report, "$where: the Lean seed ($(leanjson[:n0]), $(lean["beta"])) is not the certificate's; give lean.certificate")
    path = joinpath(root, get(lean, "certificate", rec["certificate"]))
    params = something(getk(rec, "params"), Dict())
    s, t = c.s, c.t
    if s === nothing || t === nothing
        w = window_defects(path, cert_exact(c.n0), cert_exact(c.beta), params)
        w === nothing || ((s, t) = w)
    end
    compare_defects!(report, where, (s, t), (leanjson[:s], leanjson[:t]))
    seed2 = get(leanjson, :seed2, nothing)
    if seed2 !== nothing
        w = window_defects(path, parse_exact(string(seed2[:n0])), parse_exact(string(seed2[:beta])), params)
        w === nothing && return err!(report, "$where: the certificate states no one-step entry for the second Lean seed ($(seed2[:n0]), $(seed2[:beta]))")
        compare_defects!(report, "$where: second seed", w, (seed2[:s], seed2[:t]))
    end
    return report
end

"Each certified defect must be at most the defect the Lean record assumes."
function compare_defects!(report, where, certified, assumed)
    for (k, eps, lean) in zip((:s, :t), certified, assumed)
        eps === nothing && (err!(report, "$where: the certificate states no defect $k"); continue)
        cert_exact(eps) <= parse_exact(string(lean)) ||
            err!(report, "$where: Lean assumes $k = $lean, below the certified defect")
    end
    return report
end

"""
    check_lean_interval!(report, where, rec, leanjson)

A row on a parameter interval (params `<p>_min`, `<p>_max`) needs a Lean record on the same
interval (`Record.SeededOn`, whose `interval` is `[a, b]`), and a Lean record on an interval
needs such a row.
"""
function check_lean_interval!(report, where, rec, leanjson)
    params = something(getk(rec, "params"), Dict())
    lo = [k for k in keys(params) if endswith(k, "_min")]
    interval = get(leanjson, :interval, nothing)
    if isempty(lo)
        interval === nothing || err!(report, "$where: the Lean record is on an interval but the row is a point")
        return report
    end
    length(lo) == 1 || return err!(report, "$where: more than one parameter interval")
    hi = replace(lo[1], r"_min$" => "_max")
    interval === nothing && return err!(report, "$where: the row is an interval but the Lean record is a point")
    [parse_exact(params[lo[1]]), parse_exact(get(params, hi, "0"))] == parse_exact.(string.(interval)) ||
        err!(report, "$where: the Lean interval $(collect(interval)) is not [$(lo[1]), $hi] of the row")
    return report
end

"""
    lean_covers(bound, lean, headline, leanbound, lengths) -> Bool

Whether the Lean theorem Δ_L > log(r)/β for every even L ≥ L0 implies the headline claim: the
lengths are even and periodic, and the Lean value is at least the headline value (identical
{r, beta}, or disjoint outward enclosures). Throws when the enclosures overlap.
"""
function lean_covers(bound, lean, headline, leanbound, lengths)
    something(getk(lengths, "parity"), "even") == "even" &&
        something(getk(lengths, "boundary"), "periodic") == "periodic" || return false
    same = haskey(bound, "r") && parse_exact(bound["r"]) == parse_exact(lean["r"]) &&
        parse_exact(bound["beta"]) == parse_exact(lean["beta"])
    (same || leanbound[1] >= headline[2]) && return true
    leanbound[2] < headline[1] && return false
    throw(ArgumentError("cannot decide whether the Lean bound is at least the headline bound"))
end

function check_record!(report, rec, ctx)
    (; root, families, ids, used_lean, json, summary) = ctx
    id = getk(rec, "id")
    where = id isa AbstractString ? "record $id" : "a record"
    if !(id isa AbstractString) || !occursin(ID, id)
        err!(report, "$where: id must be kebab-case")
    elseif id in ids
        err!(report, "$where: duplicate id")
    else
        push!(ids, id)
    end
    tier = getk(rec, "tier")
    tier in TIERS || return err!(report, "$where: tier $(repr(tier)) is not one of $(join(TIERS, ", "))")
    summary[tier] += 1
    fam = get(families, getk(rec, "family"), nothing)
    fam === nothing && err!(report, "$where: unknown family $(repr(getk(rec, "family")))")
    need_str!(report, where, rec, "parameter")
    slug = something(getk(rec, "slug"), getk(fam, "slug"), Some(nothing))
    params = something(getk(rec, "params"), Dict())
    if params isa AbstractDict
        foreach(((k, v),) -> check_exact!(report, "$where: params.$k", v; positive = false), params)
    else
        err!(report, "$where: params must be a mapping")
    end
    gap = getk(rec, "numerical_gap")
    gap === nothing || check_exact!(report, "$where: numerical_gap", gap)

    if tier == "target"
        for k in ("bound", "lean", "lengths")
            getk(rec, k) === nothing || err!(report, "$where: a target has no $k")
        end
        est = getk(rec, "estimate")
        est isa AbstractDict || return err!(report, "$where: a target needs an estimate")
        foreach(k -> haskey(est, k) && check_exact!(report, "$where: estimate.$k", est[k]), ("n0", "beta0", "float_bound"))
        return report
    end

    lengths = getk(rec, "lengths")
    L0 = getk(lengths, "L0")
    (L0 isa Integer && L0 > 0) || err!(report, "$where: lengths.L0 must be a positive integer")
    something(getk(lengths, "parity"), "even") in ("even", "odd", "all") || err!(report, "$where: bad lengths.parity")
    something(getk(lengths, "boundary"), "periodic") in ("periodic", "open") || err!(report, "$where: bad lengths.boundary")

    bound = getk(rec, "bound")
    headline = nothing
    if bound isa AbstractDict
        for k in keys(bound)
            k in ("r", "beta", "value", "strict", "decimal", "display") || err!(report, "$where: unknown bound field $k")
        end
        try
            headline = bound_enclosure(bound)
            check_decimal!(report, "$where: bound", getk(bound, "decimal"), headline...)
        catch e
            err!(report, "$where: bound: $(sprint(showerror, e))")
        end
    else
        err!(report, "$where: missing bound")
    end
    seed = getk(rec, "seed")
    (getk(seed, "n0") isa Integer && getk(seed, "n0") > 0) || err!(report, "$where: seed.n0 must be a positive integer")
    check_exact!(report, "$where: seed.beta0", getk(seed, "beta0"))
    check_date!(report, where, getk(rec, "date"))
    pr = getk(rec, "pr")
    (pr === nothing || (pr isa Integer && pr > 0)) || err!(report, "$where: pr must be a number or null")
    people = getk(rec, "contributors")
    if !(people isa AbstractVector) || isempty(people)
        err!(report, "$where: contributors must be a nonempty list")
    else
        for c in people
            (c isa AbstractString && (occursin(HANDLE, c) || occursin(' ', strip(c)))) ||
                err!(report, "$where: contributor $(repr(c)) is neither a GitHub handle nor a name with a space")
        end
    end
    cert = getk(rec, "certificate")
    cert === nothing ? err!(report, "$where: missing certificate") : check_path!(report, "$where: certificate", root, cert)
    check = something(getk(rec, "check"), slug === nothing ? nothing : "models/$slug/check.jl", Some(nothing))
    check === nothing ? err!(report, "$where: no checker") : check_path!(report, "$where: check", root, check)

    rowcert = cert === nothing ? nothing : check_certificate!(report, where, rec, headline, root)

    lean = getk(rec, "lean")
    if haskey(LEAN_KIND, tier) && !(lean isa AbstractDict)
        err!(report, "$where: tier $tier needs a Lean record")
    elseif lean isa AbstractDict
        name = getk(lean, "name")
        if !(name isa AbstractString) || !occursin(LEAN_NAME, name)
            err!(report, "$where: lean.name must be Leaderboard.Records.<Model>.<point>")
        else
            name in used_lean && err!(report, "$where: Lean record $name used twice")
            push!(used_lean, name)
        end
        try
            leanbound = bound_enclosure(Dict("r" => lean["r"], "beta" => lean["beta"]))
            haskey(lean, "decimal") && check_decimal!(report, "$where: lean", lean["decimal"], leanbound...)
            if headline !== nothing
                covers = lean_covers(bound, lean, headline, leanbound, lengths)
                covers && !haskey(LEAN_KIND, tier) &&
                    err!(report, "$where: the Lean theorem proves the headline bound, so the tier must be lean-unconditional or lean-conditional, not $tier")
                !covers && haskey(LEAN_KIND, tier) &&
                    err!(report, "$where: tier $tier but the Lean bound is weaker than the headline bound; use the certificate's tier")
            end
        catch e
            err!(report, "$where: lean: $(sprint(showerror, e))")
        end
        if json !== nothing && name isa AbstractString && L0 isa Integer
            kind = check_lean_json!(report, where, lean, L0, json)
            haskey(LEAN_KIND, tier) && kind isa AbstractString && kind != LEAN_KIND[tier] &&
                err!(report, "$where: tier $tier but the Lean record is $kind")
            kind == "conditional" && check_lean_seeds!(report, where, rec, lean, json[name], rowcert, root)
            haskey(json, name) && check_lean_interval!(report, where, rec, json[name])
        end
    end

    previous = nothing
    for (j, h) in enumerate(something(getk(rec, "history"), Any[]))
        hw = "$where: history[$j]"
        check_date!(report, hw, getk(h, "date"))
        need_str!(report, hw, h, "method", "lengths")
        getk(h, "tier") in TIERS || err!(report, "$hw: tier must be one of $(join(TIERS, ", "))")
        key = try
            (string(getk(h, "date")), first(parse_decimal(getk(h, "decimal"))))
        catch e
            err!(report, "$hw: $(sprint(showerror, e))")
            continue
        end
        # The page lists history in file order: newest first, larger bound first on a shared date.
        previous !== nothing && (key[1] > previous[1] || (key[1] == previous[1] && key[2] > previous[2])) &&
            err!(report, "$hw: history must be newest first, larger bound first on the same date")
        previous = key
    end
    return report
end

"""
    check_records(path; root, lean_records = nothing) -> (report, summary)

Validate `records.yml` (schema in its header). `root` is the repository root for paths;
`lean_records` is the generated `lean_records.json` (skipped when `nothing`). `summary` counts
the records of each tier.
"""
function check_records(path; root, lean_records = nothing)
    report = Report()
    summary = Dict(t => 0 for t in TIERS)
    data = try
        YAML.load_file(path)
    catch e
        return err!(report, "$path: YAML parse error: $(sprint(showerror, e))"), summary
    end
    data isa AbstractDict || return err!(report, "$path: top level must be a mapping"), summary
    for key in keys(data)
        key in ("tiers", "families", "challenges", "records") || err!(report, "unknown top-level key \"$key\"")
    end
    tiers = something(getk(data, "tiers"), Any[])
    sort([string(getk(t, "key")) for t in tiers]) == sort(collect(TIERS)) ||
        err!(report, "tiers must list exactly $(join(TIERS, ", "))")
    foreach(t -> need_str!(report, "tier $(getk(t, "key"))", t, "glyph", "label", "meaning"), tiers)

    families = Dict{String,Any}()
    for fam in something(getk(data, "families"), Any[])
        key = getk(fam, "key")
        where = "family $(repr(key))"
        isstr(key) || (err!(report, "$where: missing key"); continue)
        haskey(families, key) && err!(report, "$where: duplicate")
        families[key] = fam
        need_str!(report, where, fam, "name", "hamiltonian", "cells", "gap_source")
        slug = getk(fam, "slug")
        slug === nothing || check_path!(report, where, root, joinpath("models", slug); dir = true)
        grid = getk(fam, "grid")
        grid === nothing || need_str!(report, "$where: grid", grid, "rows", "cols")
    end

    json = nothing
    if lean_records !== nothing
        json = try
            load_lean_records(lean_records)
        catch e
            err!(report, "lean records: $(sprint(showerror, e))")
            nothing
        end
    end
    ctx = (; root, families, ids = Set{String}(), used_lean = Set{String}(), json, summary)
    records = getk(data, "records")
    records isa AbstractVector || return err!(report, "records must be a list"), summary
    foreach(rec -> check_record!(report, rec, ctx), records)

    for ch in something(getk(data, "challenges"), Any[])
        id = getk(ch, "id")
        where = "challenge $(repr(id))"
        (id isa AbstractString && occursin(ID, id) && !(id in ctx.ids)) ||
            err!(report, "$where: id must be kebab-case and distinct from other ids")
        id isa AbstractString && push!(ctx.ids, id)
        need_str!(report, where, ch, "title", "need")
        fam = getk(ch, "family")
        fam === nothing || haskey(families, fam) || err!(report, "$where: unknown family $(repr(fam))")
        getk(ch, "link") === nothing || check_link!(report, where, root, ch["link"])
        getk(ch, "check") === nothing || check_path!(report, "$where: check", root, ch["check"])
    end

    if json !== nothing
        for name in sort!(collect(setdiff(keys(json), ctx.used_lean)))
            push!(report.warnings, "Lean record $name has no entry in records.yml")
        end
        push!(report.notes, "$(length(ctx.used_lean)) Lean records matched against $(length(json)) in lean_records.json")
    end
    return report, summary
end

"""
    check_news(path; record_ids) -> Report

Validate `news.yml`: a list of `{date, text, link, prs, records}` with ISO dates, a link that is
an http(s) URL or a site path starting with `/` (optional when `prs` is given), pull request
numbers and existing record ids.
"""
function check_news(path; record_ids = Set{String}())
    report = Report()
    items = try
        YAML.load_file(path)
    catch e
        return err!(report, "$path: YAML parse error: $(sprint(showerror, e))")
    end
    items isa AbstractVector || return err!(report, "$path: expected a list of news items")
    for (i, item) in enumerate(items)
        where = "news[$i]"
        check_date!(report, where, getk(item, "date"))
        need_str!(report, where, item, "text")
        link = getk(item, "link")
        if link === nothing
            isempty(something(getk(item, "prs"), Any[])) && err!(report, "$where: needs a link or prs")
        elseif !(occursin(URL, string(link)) || occursin(r"^/[^\s\"<>]*$", string(link)))
            err!(report, "$where: link must be an http(s) URL or a site path starting with /")
        end
        for pr in something(getk(item, "prs"), Any[])
            (pr isa Integer && pr > 0) || err!(report, "$where: prs must be pull request numbers")
        end
        for id in something(getk(item, "records"), Any[])
            id in record_ids || err!(report, "$where: unknown record $(repr(id))")
        end
    end
    return report
end

const USAGE = """
usage: julia --project=tools/leaderboard tools/leaderboard/check_records.jl [options]
  --records PATH       leaderboard data (default home_page/_data/records.yml)
  --news PATH          news list (default home_page/_data/news.yml)
  --lean-records PATH  lean_records.json from `lake exe leaderboard_audit` (optional)
  --root PATH          repository root (default: two levels above this tool)
"""

"Command-line entry point; returns the process exit code."
function main(args; default_root = normpath(joinpath(@__DIR__, "..", "..", "..")))
    opts = Dict("--root" => default_root)
    i = 1
    while i <= length(args)
        a = args[i]
        a in ("-h", "--help") && (print(USAGE); return 0)
        if a in ("--records", "--news", "--lean-records", "--root") && i < length(args)
            opts[a] = args[i+1]; i += 2
        else
            print(stderr, "unknown argument $a\n", USAGE); return 2
        end
    end
    root = opts["--root"]
    records = get(opts, "--records", joinpath(root, "home_page", "_data", "records.yml"))
    report, summary = check_records(records; root, lean_records = get(opts, "--lean-records", nothing))
    haskey(opts, "--lean-records") || push!(report.notes, "no --lean-records given: Lean records not compared")
    news = get(opts, "--news", joinpath(root, "home_page", "_data", "news.yml"))
    if isfile(news) && passed(report)
        ids = Set(String[string(r["id"]) for r in YAML.load_file(records)["records"]])
        append!(report.errors, check_news(news; record_ids = ids).errors)
    end

    println("Leaderboard records: $(sum(values(summary))) in $records")
    foreach(t -> println("  ", rpad(t, 20), lpad(summary[t], 4)), TIERS)
    foreach(m -> println("note: ", m), report.notes)
    foreach(m -> println("warning: ", m), report.warnings)
    foreach(m -> println(stderr, "error: ", m), report.errors)
    passed(report) && (println("PASSED"); return 0)
    println(stderr, "FAILED: $(length(report.errors)) error(s)")
    return 1
end

end
