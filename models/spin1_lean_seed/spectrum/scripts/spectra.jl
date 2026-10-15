include("../src/common.jl")

function run_ring(n,twist)
    output=joinpath(DAT,"spectrum_$(n)_$(twist).json")
    isfile(output) && return
    records = []
    E0 = Inf
    for m in 0:n
        H,_ = hamiltonian(n,m,twist)
        # m=0 contains a ground state for these rings; verify all other sectors.
        F=low_spectrum(H;dense=n<=8,cutoff=E0+3.15,seed=100n+m)
        E0=min(E0,first(F.values))
        if m == 0 && n>8 && last(F.values) < E0+3.15
            F=low_spectrum(H;cutoff=E0+3.15,initial=64,seed=100n+m)
        end
        m == 0 && (E0=first(F.values))
        push!(records,(;m,dimension=size(H,1),nnz=nnz(H),values=F.values,
                         residual=F.residual,nmult=F.nmult,seconds=F.seconds,
                         method=F.method))
        println("spectrum n=$n twist=$twist m=$m d=$(size(H,1)) ",
                "nev=$(length(F.values)) E=$(first(F.values)) res=$(F.residual)")
        flush(stdout)
    end
    @assert E0 == first(records[1].values)
    writejson(output,merge(metadata(["src/common.jl","scripts/spectra.jl"]),
              Dict("n"=>n,"twist"=>twist,"E0"=>E0,"sectors"=>records)))
end

if abspath(PROGRAM_FILE)==@__FILE__
    ns=isempty(ARGS) ? collect(4:12) : parse.(Int,ARGS)
    for n in ns, twist in (n in (6,8,10) ? TWISTS : TWISTS[1:2])
        run_ring(n,twist)
    end
end
