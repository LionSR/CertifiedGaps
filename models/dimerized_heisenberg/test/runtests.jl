include("../check.jl")
@testset "Exact RP obstruction audit (not a gap certificate)" begin
    check_audit()
end
include("numerical_blocks.jl")
include("../reflection_conventions.jl")
@testset "Exact alternative cell conventions" begin
    conventions=ReflectionConventions.audit_conventions()
    @test length(conventions)==2
    @test all(row.reflected_form && !row.uniform_identical_form for row in conventions)
end
include("verified_independent.jl")
