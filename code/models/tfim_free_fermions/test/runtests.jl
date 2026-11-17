# Fast tests of the free-fermion TFIM certificates.
#
#     julia --project=code/julia code/models/tfim_free_fermions/test/runtests.jl
using Test
include(joinpath(@__DIR__, "..", "src", "common.jl"))
import Nemo

# Dense H_L in Float64 in the σᶻ product basis (bit j of the index is site j: 0 ↔ Z = +1),
# for the comparison with exact diagonalization.
function dense_tfim(L, g)
    H = zeros(2^L, 2^L)
    for b in 0:(2^L-1)
        z(j) = iseven(b >> mod(j, L)) ? 1 : -1
        H[b+1, b+1] = -sum(z(j) * z(j + 1) for j in 0:(L-1))
        for j in 0:(L-1)
            H[(b⊻(1<<j))+1, b+1] -= g
        end
    end
    return H
end

# Independent evaluation of the closed form in Arb balls (Nemo), straight from the four
# products without the rescaling of src/FreeFermions.jl.
function arb_logZ(L, g, β; bits = 1024)
    RR = Nemo.ArbField(bits)
    # Rationals from their integers: Nemo converts a Julia Rational through a BigFloat at the
    # global precision, which is inexact for 11/10.
    q(x) = RR(numerator(x)) / RR(denominator(x))
    G, B = q(g), q(β)
    e(m, s) = 2 * sqrt(1 + G^2 - 2G * cos(2 * Nemo.const_pi(RR) * (m + q(s)) / L))
    prods(s) = (prod(2cosh(B * e(m, s) / 2) for m in 0:(L-1)),
        prod(2sinh(B * e(m, s) / 2) for m in 0:(L-1)))
    (cn, sn), (cr, sr) = prods(1 // 2), prods(0)
    return log((cn + sn + cr - sr) / 2)
end

@testset "closed form vs exact diagonalization" begin
    for L in 4:2:10, g in (11 // 10, 5 // 4, 2 // 1, 3 // 1)
        ev = eigvals(Symmetric(dense_tfim(L, Float64(g))))
        # Ground state: the filled Neveu–Schwarz sea, E₀ = −½ Σ_NS ε_k.
        @test ev[1] ≈ -sum(2sqrt(1 + g^2 - 2g * cospi((2m + 1) / L)) for m in 0:(L-1)) / 2
        for β in (1 // 3, 2 // 1, 7 // 1)
            ed = log(sum(exp.(-Float64(β) .* (ev .- ev[1])))) - Float64(β) * ev[1]
            x = logZ(L, g, big(β))
            @test isapprox(Float64(mid(x)), ed; rtol = 1e-12)
        end
    end
end

@testset "closed form vs Arb referee" begin
    for (L, g, β) in ((16, 2 // 1, 3 // 1), (64, 5 // 4, 21 // 2), (128, 11 // 10, 40 // 1))
        x, ref = logZ(L, g, big(β)), arb_logZ(L, g, β)
        RR = Nemo.ArbField(1024)
        lo, hi = (RR(numerator(exact(e))) / RR(denominator(exact(e))) for e in (inf(x), sup(x)))
        @test lo <= ref && ref <= hi
        @test diam(x) < big(10)^-100
    end
end

@testset "domain, RP record, and BigFloat precision traps" begin
    x = setprecision(() -> big(1) / 3, BigFloat, 512)
    @test BigFloat(exact(x); precision = 512) == x
    @test_throws ArgumentError logZ(8, 1 // 1, big(1) // 1)
    @test_throws ArgumentError logZ(8, 1 // 2, big(1) // 1)
    @test rp_record(2 // 1)["C"] == [["1", "0"], ["0", "-1"]]
    # seed_enclosures must give enclosures of width > 0 at PRECISION bits (see common.jl).
    s = seed_enclosures(5 // 4, 8, 2 // 1)
    @test 0 < s.eps_s < 1 && 0 < s.eps_t < 1
    ln, l2, l22 = (Float64(mid(x)) for x in s.logZ)
    @test isapprox(Float64(s.eps_s), -expm1(l2 - 2ln); rtol = 1e-6)
end

@testset "stored certificates" begin
    for p in POINTS
        data = read_result(certificate_path(p.g), Dict{String,Any})
        b, c = check_point(data, p.g, p.n0, p.beta0)
        @test 0 < c < b < 2(p.g - 1)
    end
end

@testset "tampering is rejected" begin
    p = POINTS[1]
    data = read_result(certificate_path(p.g), Dict{String,Any})
    # A smaller εₛ changes the rerun bound, so the stored bound no longer matches.
    bad = deepcopy(data)
    bad["seed"]["eps_s"] = string(parse_rational(data["seed"]["eps_s"]) / 2)
    @test_throws ErrorException check_point(bad, p.g, p.n0, p.beta0)
    # With the bound also forged to match, only the full replay of log Z catches it.
    ib, chain = bounds(parse_rational(bad["seed"]["eps_s"]),
        parse_rational(bad["seed"]["eps_t"]), p.n0, p.beta0)
    bad["seed"]["bound_exact_binary"] = string(RationalBig(ib.bound))
    bad["chain"] = chain
    @test check_point(bad, p.g, p.n0, p.beta0) isa Tuple
    @test_throws ErrorException check_point(bad, p.g, p.n0, p.beta0; full = true)
    @test_throws ErrorException check_point(data, 5 // 4, p.n0, p.beta0)
end
