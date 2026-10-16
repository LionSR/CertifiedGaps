# Certificate verification, shared by `check.jl` (the entry point) and `test/runtests.jl`.
# Included by `LongitudinalCertificate`.

"""
    check_result(path; partitions=true, full=false)

Verify one per-point certificate and return `true`; throw on any mismatch. With
`partitions=false` the (slower) interval partition sums are skipped. With `full=true` every
block is rebuilt with its rigorous enclosure, and each fresh interval must lie, in exact
rationals, inside the stored one, so the stored radii are re-proved rather than trusted.
"""
function check_result(path; partitions=true, full=false)
    d=read_result(path,Dict{String,Any}); g=parse_rational(d["g"]); h=parse_rational(d["h"])
    @check rp_check(g,h)
    haskey(d,"certificate") || error("$(basename(path)): no certificate stored")
    c=d["certificate"]; n=c["n"]; @check iseven(n)
    for (key,L) in (("small",n),("large",2n))
        @check sum(b["dimension"] for b in c[key])==2^L
        @check all(length(b["values"])==b["dimension"] for b in c[key])
        @check [b["dimension"] for b in c[key]]==
            filter(!iszero,sector_dimensions(orbits(L)))
    end
    if full
        for (key,L) in (("small",n),("large",2n))
            fresh=ring_spectrum(L,g,h;rigorous=true)
            @check length(fresh.blocks)==length(c[key])
            for (new,old) in zip(fresh.blocks,c[key])
                @check new["k"]==old["k"] && new["dimension"]==old["dimension"]
                # The theorem proves radius/2; stored intervals include a factor-two margin.
                fresh_radius=parse_rational(new["radius"])/2
                stored_radius=parse_rational(old["radius"])
                @check all(abs(parse_rational(a)-parse_rational(b))+fresh_radius<=stored_radius
                    for (a,b) in zip(new["values"],old["values"]))
            end
        end
    end
    b=parse_rational(c["beta"])
    es=parse_rational(c["epsilon_s"]); et=parse_rational(c["epsilon_t"])
    if partitions
        seeds=seed_enclosure(c["small"],c["large"],b)
        @check seeds.es==es && seeds.et==et
    end
    bound=certified_gap_bound(es,et,b).bound
    @check bound!==nothing && bound>0
    @check bound==c["bound"]
    @check exact(bound)==parse_rational(c["bound_exact_binary"])
    true
end

"Check that the stored certificates are exactly the 30 grid points `GS × HS`."
function check_grid(paths)
    @check length(paths)==length(GS)*length(HS)
    points=Set(let d=read_result(p,Dict{String,Any})
        (parse_rational(d["g"]),parse_rational(d["h"])) end for p in paths)
    @check points==Set((g,h) for g in GS for h in HS)
    true
end
