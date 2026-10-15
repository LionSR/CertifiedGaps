# Consistency checks of the stored numerical data. This is NOT a certificate checker: the
# spectra are Float64 (dense LAPACK / ARPACK) and nothing here proves a spectral count.
#
#   julia --project=. check.jl      # about 10 s; writes results/reference_trace_checks.csv and
#                                   # the (git-ignored) results/dense_validation.csv
isdefined(@__MODULE__, :DAT) || include("src/common.jl")
using Test

"""
    check_hamiltonian()

Recompute the Hamiltonians and compare with the stored spectra: basis sizes, traces, forced
ARPACK against the stored dense spectra through n = 8, and m <-> -m spectral equality.
"""
function check_hamiltonian()
    rows=[]
    @testset "Hamiltonian and multiplicities" begin
        for n in 4:12
            @test sum(length(basis(n,m)) for m in -n:n)==3^n
        end
        for n in 4:8, twist in (n in (6,8) ? TWISTS : TWISTS[1:2])
            ring=JSON3.read(read(joinpath(DAT,"spectrum_$(n)_$(twist).json"),String))
            for s in ring.sectors
                H,_=hamiltonian(n,Int(s.m),twist)
                @test tr(H)==sum(Float64.(s.values)) ||
                      abs(tr(H)-sum(Float64.(s.values)))<1e-9
                size(H,1)>40 || continue
                F=low_spectrum(H;initial=min(24,size(H,1)-2),seed=271+n,force_sparse=true)
                err=maximum(abs.(F.values.-Float64.(s.values[1:length(F.values)])))
                @test err<1e-8
                push!(rows,(;n,twist,m=Int(s.m),nev=length(F.values),error=err,
                            residual=F.residual))
            end
        end
        for twist in TWISTS
            Hp,_=hamiltonian(6,1,twist); Hm,_=hamiltonian(6,-1,twist)
            @test eigvals(Hermitian(Matrix(Hp))) ≈ eigvals(Hermitian(Matrix(Hm)))
        end
        @test first(eigvals(Hermitian(Matrix(first(hamiltonian(4,0,"I")))))) ≈ -6
    end
    writecsv(joinpath(DAT,"dense_validation.csv"),rows)
end

"""
    check_data()

Check the stored spectra (dimensions, residuals, guard levels, trace and second moment through
n = 8), the shifted traces against the centres of `models/spin1/results/source{1,2}.json`, and
that every adaptive filter schedule in summary.json met its target.
"""
function check_data()
    analysis=JSON3.read(read(joinpath(DAT,"analysis.json"),String))
    summary=JSON3.read(read(joinpath(DAT,"summary.json"),String))
    rows=[]
    @testset "Complete numerical experiment" begin
        @test length(analysis.traces)==42
        for n in 4:12, twist in (n in (6,8,10) ? TWISTS : TWISTS[1:2])
            f=JSON3.read(read(joinpath(DAT,"spectrum_$(n)_$(twist).json"),String))
            @test sum(s.dimension*(s.m==0 ? 1 : 2) for s in f.sectors)==3^n
            @test maximum(s.residual for s in f.sectors)<1e-8
            for s in f.sectors
                @test length(s.values)==s.dimension || last(s.values)>f.E0+3.15
            end
            if n<=8
                @test abs(sum(sum(s.values)*(s.m==0 ? 1 : 2) for s in f.sectors))<1e-8
                @test sum(sum(abs2,s.values)*(s.m==0 ? 1 : 2) for s in f.sectors) ≈
                      4n*3^(n-1)
            end
        end
        for source in (1,2)
            path=joinpath(ROOT,"../../spin1/results/source$source.json")
            reference=JSON3.read(read(path,String))
            for c in reference.certificates
                twist=TWISTS[Int(c.twist)+1]; beta=source==1 ? 12.25 : 10.5
                t=only(filter(t->t.n==c.n && t.twist==twist && t.beta==beta,analysis.traces))
                terms=split(String(c.center),"//")
                center=parse(Float64,terms[1])/parse(Float64,terms[2])
                @test abs(t.Z_shifted-center)<3e-5
                push!(rows,(;source,n=Int(c.n),twist,center,
                            measured=t.Z_shifted,error=t.Z_shifted-center))
            end
            expected=source==1 ? 44 : 24
            for bits in (64,128),method in ("linear","quadratic")
                r=filter(r->r.source==source && r.bits==bits && r.degree==0 &&
                            r.method==method,summary.runs)
                @test length(r)==expected
                @test all(x->x.success,r)
                # A small residual alone cannot exclude a missed lower state;
                # independently compare the Ritz values with the ED/Lanczos list.
                @test maximum(x.eigerror for x in r)<1e-5
            end
        end
    end
    writecsv(joinpath(DAT,"reference_trace_checks.csv"),rows)
end

if abspath(PROGRAM_FILE)==@__FILE__
    check_hamiltonian()
    check_data()
end
