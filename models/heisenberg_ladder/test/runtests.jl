using Test, LinearAlgebra, IntervalArithmetic
include("../src/Suite.jl")
using .LadderSuite
const LS = LadderSuite

@testset "Exact RP and original-to-rotated identities" begin
    for J in LS.COUPLINGS
        @test LS.exact_rp(J)
        model,rp = LS.rp_model(J)
        @test LS.verify_rp(model,rp)
        @test LS.decode_model(LS.encode_model(model)).bond == model.bond
        bad = deepcopy(LS.encode_rp(rp))
        bad["factors"][1]["weight"] = "2//1"
        @test !LS.verify_rp(model,LS.decode_rp(bad))
    end
end

@testset "Sector reconstruction and complete small spectra" begin
    for L in (2,3,4)
        J = 7//2
        H = zeros(4^L,4^L)
        for s in 0:4^L-1, (t,h) in LS.transitions(s,L,J)
            H[t+1,s+1] = Float64(h)
        end
        # Independent shared-core graph builder checks the bit transition convention.
        x = [0.0 0.5; 0.5 0.0]
        a = [0.0 0.5; -0.5 0.0]
        z = [0.5 0; 0 -0.5]
        bond = kron(x,x)-kron(a,a)+kron(z,z)
        bonds = [(2j+1,2j+2,Float64(J)*bond) for j in 0:L-1]
        append!(bonds,[(2j+b,2mod(j+1,L)+b,bond) for j in 0:L-1 for b in 1:2])
        chain = LS.Chain(2,_ -> bonds)
        @test H == Matrix(LS.sector_hamiltonian(chain,2L,collect(0:4^L-1)))
        full,rad,_ = LS.verified_eigen(H,0.0)
        reduced = LS.full_spectrum(L,J; verified=true)
        @test all(abs.(full-reduced.energies) .<= rad .+ reduced.radii)
        @test sum(d*m for (d,m) in values(LS.sector_dimensions(L))) == 4^L
        # L=2 has doubled periodic leg bonds, so its second trace moment differs.
        expected = big(4)^L*(3//16)*(L*J^2+(L==2 ? 8 : 2L))
        total = sum((interval(e)+interval(-r,r))^2
            for (e,r) in zip(reduced.energies,reduced.radii))
        @test LS.Q(inf(total)) <= expected <= LS.Q(sup(total))
    end
end

@testset "Recurrence rejection is a positive proof" begin
    @test LS.rejection(big(1)//2,big(0)//1) !== nothing
    @test LS.rejection(big(0)//1,big(1)//2) !== nothing
    @test LS.rejection(big(1)//100,big(1)//100) === nothing
    @test_throws ErrorException LS.rejection(big(-1)//100,big(1)//100)
    # Exact rational inequalities used in the arbitrary-size perturbation proof.
    @test 1+1+1//2 > 2
    @test 4/(1+9//4)^2 == 64//169 < 1//2
end

@testset "Certified coupling continuation" begin
    seed = Dict("n0"=>"4//1", "beta0"=>"2//1",
        "eps_s_upper"=>"3508321//100000000", "eps_t_upper"=>"123121//5000000")
    c = LS.continuation(seed,big(1)//10000)
    @test c !== nothing
    @test LS.q(c["eps_s_upper"]) > LS.q(seed["eps_s_upper"])
    @test LS.q(c["eps_t_upper"]) > LS.q(seed["eps_t_upper"])
    @test LS.q(c["gap_lower_bound"]) > 0
end

@testset "Continuous-temperature negative certificate and coverage" begin
    specs = Dict(L=>LS.full_spectrum(L,2//1; verified=true) for L in (2,4))
    record = LS.exclude_temperatures(specs,2)
    @test LS.check_exclusion(specs,record)
    bad = deepcopy(record)
    bad["windows"][1]["beta_from"] = "1//100"
    @test_throws ErrorException LS.check_exclusion(specs,bad)
    bad = deepcopy(record)
    bad["tail"]["eps_s_lower"] = "0//1"
    @test_throws ErrorException LS.check_exclusion(specs,bad)
end
