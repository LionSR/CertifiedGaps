include("certify.jl")
using .LongitudinalCertificate, CertifiedGaps, JSON3, CertifiedGaps.Results
const C=LongitudinalCertificate
function check_result(path; partitions=true)
    d=C.read_point(path); g=C.rat(d["g"]); h=C.rat(d["h"])
    @assert C.rp_check(g,h)
    haskey(d,"certificate") || return false
    c=d["certificate"]; n=c["n"]; @assert iseven(n)
    for (key,L) in (("small",n),("large",2n))
        @assert sum(b["dimension"] for b in c[key])==2^L
        @assert all(length(b["values"])==b["dimension"] for b in c[key])
        @assert [b["dimension"] for b in c[key]]==filter(!iszero,C.sector_dimensions(C.orbits(L)))
    end
    if get(ENV,"DEEP_CHECK","0")=="1"
        for (key,L) in (("small",n),("large",2n))
            fresh=C.ring_spectrum(L,g,h;rigorous=true)
            @assert length(fresh.blocks)==length(c[key])
            for (new,old) in zip(fresh.blocks,c[key])
                @assert new["k"]==old["k"] && new["dimension"]==old["dimension"]
                # The theorem proves radius/2; stored intervals include a factor-two margin.
                fresh_radius=C.rat(new["radius"])/2
                stored_radius=C.rat(old["radius"])
                @assert all(abs(C.rat(a)-C.rat(b))+fresh_radius<=stored_radius
                    for (a,b) in zip(new["values"],old["values"]))
            end
        end
    end
    b=C.rat(c["beta"]); es=C.rat(c["epsilon_s"]); et=C.rat(c["epsilon_t"])
    if partitions
        seeds=C.seed_enclosure(c["small"],c["large"],b)
        @assert seeds.es==es && seeds.et==et
    end
    bound=certified_gap_bound(es,et,b).bound
    @assert bound==c["bound"]
    @assert C.exact(bound)==C.rat(c["bound_exact_binary"])
    true
end
if abspath(PROGRAM_FILE)==@__FILE__
    for path in C.result_paths()
        println(basename(path),": ",check_result(path))
    end
end
