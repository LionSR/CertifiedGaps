include(joinpath(@__DIR__, "setup.jl"))
using Test, LinearAlgebra, Random
const MN = CM.Norms
const F = MN.F

@testset "Exact dyadic tokens and legacy rational witnesses" begin
    values = F[0, 1, -1, 1//3, -7//15, 7//(big(1)<<64), -13//(big(1)<<4096)]
    for x in values
        @test MN.parse_fraction(MN.encode_fraction(x)) == x
        @test MN.parse_fraction(string(x)) == x
        @test length(MN.encode_fraction(x)) <= length(string(x))
    end
    @test MN.encode_fraction(values[end-1]) == "7//2^64"
    @test MN.parse_fraction("3//2^0") == 3
    @test MN.parse_fraction("3//-7") == -3//7
    @test_throws AssertionError MN.parse_fraction("1//2^-1")
    @test_throws AssertionError MN.parse_fraction("1//0")
    @test_throws AssertionError MN.parse_fraction("1//2//3")
    @test_throws ArgumentError MN.parse_fraction("1//2^x")
    @test_throws ArgumentError MN.parse_fraction("1//2^3.0")
    core = reshape(vcat(values, F[0]), 1, 2, 2, 2)
    legacy = [
        Dict(
            "dimensions"=>[1, 2, 2, 2],
            "entries"=>[[i, string(core[i])] for i in eachindex(core) if !iszero(core[i])],
        ),
    ]
    @test MN.decode_cores(legacy) == [core]
    @test MN.decode_cores(MN.encode_cores([core])) == [core]
    gauges = [reshape(F[1, 0, values[end-1], 1], 2, 2)]
    @test MN.decode_gauges(MN.encode_gauges(gauges)) == gauges
end

function explicit_prefixes(cores)
    operators = [ones(F, 1, 1)]
    prefixes = [operators]
    for core in cores
        l, d, _, r = size(core)
        next = [zeros(F, d*size(first(operators), 1), d*size(first(operators), 2)) for _ in 1:r]
        for a in 1:r, b in 1:l
            next[a] += kron(operators[b], core[b, :, :, a])
        end
        operators = next
        push!(prefixes, operators)
    end
    return prefixes
end

@testset "Interval local-map norm witnesses" begin
    W = F[3//4 0; 0 1//4]
    @test MN.verify_local_bound(W, F(751, 1000))
    @test !MN.verify_local_bound(W, F(749, 1000))
    @test !MN.verify_local_bound(W, F(-1))
    @test MN.verify_local_bound(zeros(F, 3, 2), F(0))
    @test !MN.verify_local_bound(W, F(0))
    @test F(3, 4) < MN.local_norm_upper(W) < F(751, 1000)
    @test MN.local_norm_upper(zeros(F, 3, 2)) == 0
    tiny = F(1, big(10)^100)
    @test tiny*F(3, 4) < MN.local_norm_upper(tiny*W) < tiny*F(751, 1000)
    for W in (F[1 2 3; 4 5 6], F[1 2; 3 4; 5 6])
        alpha = MN.local_norm_upper(W)
        @test Float64(alpha) >= opnorm(Float64.(W))
        @test !MN.verify_local_bound(W, alpha*999//1000)
    end
    rng = MersenneTwister(992)
    for dims in ([1, 2, 1], [1, 3, 2, 1], [1, 2, 3, 2, 1])
        cores = [F.(rand(rng, -3:3, dims[j], 2, 2, dims[j+1])) ./ 4 for j in 1:(length(dims)-1)]
        x = randn(rng, sum(dims[2:(end-1)]))
        gauges = MN.rational_gauges(cores, x)
        transformed = MN.gauge_transform(cores, gauges)
        exact = only(last(explicit_prefixes(cores)))
        @test only(last(explicit_prefixes(transformed))) == exact
        certificate = MN.local_product_witness(cores, gauges)
        @test MN.decode_cores(MN.encode_cores(cores)) == cores
        @test MN.verify_local_product(cores, gauges, certificate.bounds)
        corrupted = copy(certificate.bounds)
        corrupted[1] = 0
        @test !MN.verify_local_product(cores, gauges, corrupted)
        @test Float64(certificate.upper) >= opnorm(Float64.(exact))
        for (c, alpha) in zip(transformed, certificate.bounds)
            @test MN.verify_local_bound(MN.local_operator(c), alpha)
        end
    end
end

@testset "Exact MPO norm witnesses" begin
    for x in (F(0), F(1), F(2), F(1, 3), F(123456789, 1234))
        @test MN.upper_sqrt(x)^2 >= x
    end
    @test MN.upper_sqrt(F(1)) == 1
    identity = reshape(Matrix{F}(I, 5, 5), 1, 5, 5, 1)
    long = MN.witness(fill(identity, 128))
    @test long.upper == 1 && all(==(1), long.cuts)
    rng = MersenneTwister(732)
    for dims in ([1, 2, 1], [1, 3, 2, 1], [1, 2, 3, 2, 1])
        cores = [F.(rand(rng, -3:3, dims[j], 2, 2, dims[j+1])) ./ 4 for j in 1:(length(dims)-1)]
        certificate = MN.witness(cores)
        prefixes = explicit_prefixes(cores)
        for j in eachindex(prefixes), a in eachindex(prefixes[j])
            @test Float64(certificate.left[j][a])+1e-12 >= opnorm(Float64.(prefixes[j][a]))
        end
        actual = opnorm(Float64.(only(last(prefixes))))
        @test all(Float64(c)+1e-12 >= actual for c in certificate.cuts)
        quantized = MN.quantize(Float64.(c) for c in cores)
        @test quantized == cores
    end
end
