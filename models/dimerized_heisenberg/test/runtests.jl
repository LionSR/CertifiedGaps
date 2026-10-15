# Fast tests (about 20 s): the full checker (check.jl) plus independent validations of the
# sector code, the reflection conventions and the enclosures on small rings.
include(joinpath(@__DIR__, "..", "check.jl"))
include(joinpath(@__DIR__, "..", "src", "reflection_conventions.jl"))
@testset "Exact RP audit (the RP hypothesis; seeds are checked below)" begin
    check_audit()
end
@testset "Exact alternative cell conventions" begin
    conventions=ReflectionConventions.audit_conventions()
    @test length(conventions)==2
    @test all(row.reflected_form && !row.uniform_identical_form for row in conventions)
end
include("numerical_blocks.jl")
include("verified_independent.jl")
@testset "Verified saved spectral and purity certificates" begin
    check_verified_results()
end
