# Proposals may come from floating eigensolves. Every accepted finite inequality below
# is recomputed with Rational{BigInt}, or bounded Int128 arithmetic. Not a Lean checker.
include("certify.jl")

function check(path=joinpath(@__DIR__,"certificate.json"))
    data=read_exact_json(path)
    for p in data["parameters"]
        n=parse(Int,p["n"]);u=rat(p["u"]);v=rat(p["v"]);c=rat(p["c"])
        localrows=check_local(u,v,c)
        @assert localrows==p["local_psd"]
        for old in p["blocks"]
            m=parse(Int,old["Sz"])
            X=A.block(n÷2,m)
            mat=QQ.(X.H)+u*X.L+v*X.Q
            Xneg=A.block_matrices(n÷2,-m)
            @assert QQ.(Xneg.H)+u*Xneg.L+v*Xneg.Q==mat[end:-1:1,end:-1:1]
            cert=exact_spectrum_enclosure(mat)
            @assert all(rat.(old["lower"]).<=cert.lower)
            @assert all(cert.upper.<=rat.(old["upper"]))
        end
        @assert exact_counts(p["blocks"],n,c)==p["sector_anchors"]
        if haskey(p,"binned_tail")
            @assert exact_binned_tail(p["blocks"],n,c)==p["binned_tail"]
        end
        println("PASS exact n=$n $(p["kind"]): local PSD, twists, half spectra, counts")
    end
    println("PASS all finite checks. No Lean kernel certification is claimed.")
end
check()
