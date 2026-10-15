# Independent exact replay in Julia is a diagnostic; Lean checks are the trust endpoint.
using JSON, SHA
include("certify.jl")

function verify_certificate(c)
    @assert exactint(c["scale"])==SCALE
    @assert c["width"]==WIDTH && exactint(c["bias"])==BIAS
    for p in c["points"]
        n=p["n"]; d=1<<n
        V=map(p["packed"]) do x
            packed=exactint(x)
            [((packed>>(WIDTH*s))&((big(1)<<WIDTH)-1))-BIAS for s in 0:d-1]
        end
        @assert exactstr.(pack.(V))==p["packed"] # rejects unconsumed high digits
        lam=exactint.(p["lambda"])
        G=[sum(V[i].*V[j]) for i in 1:2,j in 1:2]
        @assert exactstr.(vec(G))==p["gram"]
        R=[SCALE.*apply_h(n,V[i]).-lam[i].*V[i] for i in 1:2]
        rs=[sum(r.^2) for r in R]
        @assert exactstr.(rs)==p["residual_sq"]
        @assert maximum(abs.(G-SCALE^2*I))*big(10)^12 <= SCALE^2
        @assert sum(rs)*big(10)^18<=SCALE^4
        anchor=exactint(p["anchor_20H"])
        @assert anchor == -26n+82
        @assert count(s->diagonal(n,s)<anchor+5n,0:d-1)==2
        @assert lam[1]+SCALE<lam[2] && lam[2]+SCALE<anchor*SCALE
    end
    true
end

if abspath(PROGRAM_FILE)==(@__FILE__)
    @assert verify_certificate(JSON.parsefile(joinpath(OUT,"certificate.json")))
    println("Exact Julia replay passed; this does not certify RPRing spectra.")
end
