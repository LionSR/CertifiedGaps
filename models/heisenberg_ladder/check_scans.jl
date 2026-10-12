# Validate the usable numerical evidence retained after compaction; not a certification.
include("ladder.jl")
using JSON3, Test
function check_scans()
    start=time(); parameters=Any[]
    for tag in ("4p0","3p5","3p0","2p5","2p0","1p5")
        scan=JSON3.read(read(joinpath(@__DIR__,"results","scan_J$(tag).json"),String))
        logs=Dict{Int,Dict{Rational{Int},Float64}}()
        for L in (2,3,4,5,6,8,10)
            s=JSON3.read(read(joinpath(@__DIR__,"results","spectrum_J$(tag)_L$(L).json"),String))
            @test sum(r.dimension*r.multiplicity for r in s.sectors)==4^L==s.dimension
            labels=[(r.m,r.k,r.leg_parity,r.spin_parity) for r in s.sectors]
            @test length(unique(labels))==length(labels)
            @test length(s.logZ_grid)==128
            logs[L]=Dict(parse(Rational{Int},r.beta)=>Float64(r.logZ) for r in s.logZ_grid)
        end
        defects=(n,b)->begin
            beta=Rational{Int}(b)
            (-expm1(logs[2n][beta]-2logs[n][beta]),
             -expm1(logs[2n][2beta]-2logs[2n][beta]))
        end
        rows=scan_seeds(defects,2:5,0.25:0.25:16)
        @test length(rows)==length(scan.grid)==256
        for (r,s) in zip(rows,scan.grid)
            @test r.n₀==s.n0 && r.β₀==s.beta
            @test isapprox(r.εₛ,s.eps_s;atol=5e-10,rtol=0)
            @test isapprox(r.εₜ,s.eps_t;atol=5e-10,rtol=0)
        end
        b=best_seed(filter(r->iseven(r.n₀),rows))
        @test (b===nothing)==(scan.best===nothing)
        if b!==nothing
            @test b.n₀==scan.best.n0 && b.β₀==scan.best.beta
            @test isapprox(b.bound,scan.best.bound;atol=5e-9,rtol=0)
        end
        push!(parameters,Dict("J"=>scan.J,"largest_dimension"=>4^10,
            "all_sector_sums_passed"=>true,"compacted_grid_reproduced"=>true))
    end
    write_result(joinpath(@__DIR__,"results/scan_validation.json"),Dict(
        "metadata"=>provenance(Dict("scope"=>"numerical artifact consistency only")),
        "parameters"=>parameters,"seconds"=>time()-start))
    println("All six compacted numerical grids and exact sector counts replay correctly.")
end
check_scans()
