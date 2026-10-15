# Checker for the 30 longitudinal-field TFIM certificates in results/g*_h*.json.gz.
#
# Fast mode (default, ~20 s): exact RP identity, momentum-sector dimensions, complete spectra,
# the interval partition sums and rational seeds recomputed from the stored eigenvalue
# intervals, and the recurrence bound. It trusts the stored eigenvalue radii.
# DEEP_CHECK=1 (~6 min): also rebuilds every block and its rigorous enclosure with the
# current code and checks, in exact rationals, that each fresh interval lies inside the stored
# one, so the radii are re-proved rather than trusted. Any failure throws.
include(joinpath(@__DIR__, "src", "LongitudinalCertificate.jl"))
using .LongitudinalCertificate, CertifiedGaps, JSON3, CertifiedGaps.Results
using .LongitudinalCertificate: @check
const C=LongitudinalCertificate

"""
    check_result(path; partitions=true)

Verify one per-point certificate and return `true`; throw on any mismatch. With
`partitions=false` the (slower) interval partition sums are skipped.
"""
function check_result(path; partitions=true)
    d=C.read_point(path); g=C.rat(d["g"]); h=C.rat(d["h"])
    @check C.rp_check(g,h)
    haskey(d,"certificate") || error("$(basename(path)): no certificate stored")
    c=d["certificate"]; n=c["n"]; @check iseven(n)
    for (key,L) in (("small",n),("large",2n))
        @check sum(b["dimension"] for b in c[key])==2^L
        @check all(length(b["values"])==b["dimension"] for b in c[key])
        @check [b["dimension"] for b in c[key]]==
            filter(!iszero,C.sector_dimensions(C.orbits(L)))
    end
    if get(ENV,"DEEP_CHECK","0")=="1"
        for (key,L) in (("small",n),("large",2n))
            fresh=C.ring_spectrum(L,g,h;rigorous=true)
            @check length(fresh.blocks)==length(c[key])
            for (new,old) in zip(fresh.blocks,c[key])
                @check new["k"]==old["k"] && new["dimension"]==old["dimension"]
                # The theorem proves radius/2; stored intervals include a factor-two margin.
                fresh_radius=C.rat(new["radius"])/2
                stored_radius=C.rat(old["radius"])
                @check all(abs(C.rat(a)-C.rat(b))+fresh_radius<=stored_radius
                    for (a,b) in zip(new["values"],old["values"]))
            end
        end
    end
    b=C.rat(c["beta"]); es=C.rat(c["epsilon_s"]); et=C.rat(c["epsilon_t"])
    if partitions
        seeds=C.seed_enclosure(c["small"],c["large"],b)
        @check seeds.es==es && seeds.et==et
    end
    bound=certified_gap_bound(es,et,b).bound
    @check bound!==nothing && bound>0
    @check bound==c["bound"]
    @check C.exact(bound)==C.rat(c["bound_exact_binary"])
    true
end

"Check that the stored certificates are exactly the 30 grid points `GS × HS`."
function check_grid(paths)
    @check length(paths)==length(C.GS)*length(C.HS)
    points=Set(let d=C.read_point(p); (C.rat(d["g"]),C.rat(d["h"])) end for p in paths)
    @check points==Set((g,h) for g in C.GS for h in C.HS)
    true
end

if abspath(PROGRAM_FILE)==@__FILE__
    paths=C.result_paths()
    check_grid(paths)
    for path in paths
        println(basename(path),": ",check_result(path))
    end
    println("All $(length(paths)) certificates verified",
        get(ENV,"DEEP_CHECK","0")=="1" ? " (deep replay)." : ".")
end
