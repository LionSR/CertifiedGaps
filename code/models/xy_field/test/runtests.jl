# Fast tests of the XY-in-a-field certificates.
#
#     julia --project=code/julia code/models/xy_field/test/runtests.jl
using Test
include(joinpath(@__DIR__, "..", "src", "common.jl"))
import Nemo

# Dense H_L in the σᶻ product basis (bit j of the index is site j: 0 ↔ Z = +1). X_j X_k flips
# both spins; Y_j Y_k does too, with the factor (i z_j)(i z_k) = −z_j z_k.
function dense_xy(L, Jx, Jy, h)
    H = zeros(2^L, 2^L)
    for b in 0:(2^L-1)
        z(j) = iseven(b >> mod(j, L)) ? 1 : -1
        H[b+1, b+1] = -h * sum(z(j) for j in 0:(L-1))
        for j in 0:(L-1)
            c = b ⊻ (1 << j) ⊻ (1 << mod(j + 1, L))
            H[c+1, b+1] += Jx - Jy * z(j) * z(j + 1)
        end
    end
    return H
end

# The same Hamiltonian from Kronecker products of Pauli matrices, as a check of dense_xy.
function kron_xy(L, Jx, Jy, h)
    X, Y, Z, I2 = ComplexF64[0 1; 1 0], ComplexF64[0 -im; im 0], ComplexF64[1 0; 0 -1],
        Matrix{ComplexF64}(I, 2, 2)
    # Site 0 is the least significant bit, i.e. the last Kronecker factor.
    op(ops) = reduce(kron, [get(ops, j, I2) for j in (L-1):-1:0])
    return sum(Jx * op(Dict(j => X, mod(j + 1, L) => X)) + Jy * op(Dict(j => Y, mod(j + 1, L) => Y))
               - h * op(Dict(j => Z)) for j in 0:(L-1))
end

# Independent evaluation of the closed form in Arb balls (Nemo), straight from the four
# products. Rationals are built from their integers: Nemo converts a Julia Rational through a
# BigFloat at the global precision.
function arb_logZ(L, Jx, Jy, h, β; bits = 1024)
    RR = Nemo.ArbField(bits)
    q(x) = RR(numerator(x)) / RR(denominator(x))
    A, B, H, b = q(Jx + Jy), q(Jx - Jy), q(h), q(β)
    function e(m, s)
        k = 2 * Nemo.const_pi(RR) * (m + q(s)) / L
        return 2 * sqrt((H - A * cos(k))^2 + B^2 * sin(k)^2)
    end
    prods(s) = (prod(2cosh(b * e(m, s) / 2) for m in 0:(L-1)),
        prod(2sinh(b * e(m, s) / 2) for m in 0:(L-1)))
    (cn, sn), (cr, sr) = prods(1 // 2), prods(0 // 1)
    return log((cn + sn + cr - sr) / 2)
end

const PARAMS = [(-1 // 1, 1 // 4, 1 // 1), (-1 // 1, 1 // 2, 1 // 1), (-1 // 1, 1 // 1, 1 // 1),
    (-1 // 1, 1 // 2, 2 // 1), (1 // 2, -1 // 1, 3 // 2)]

@testset "dense Hamiltonian" begin
    for (Jx, Jy, h) in PARAMS
        @test dense_xy(4, Float64.((Jx, Jy, h))...) ≈ kron_xy(4, Float64.((Jx, Jy, h))...)
    end
end

@testset "closed form vs exact diagonalization" begin
    for L in 4:2:10, (Jx, Jy, h) in PARAMS
        ev = eigvals(Symmetric(dense_xy(L, Float64.((Jx, Jy, h))...)))
        # Ground state: the filled Neveu–Schwarz sea.
        a, b = Jx + Jy, Jx - Jy
        @test ev[1] ≈ -sum(sqrt((h - a * cospi((2m + 1) / L))^2 + b^2 * sinpi((2m + 1) / L)^2)
                           for m in 0:(L-1))
        for β in (1 // 3, 2 // 1, 7 // 1)
            ed = log(sum(exp.(-Float64(β) .* (ev .- ev[1])))) - Float64(β) * ev[1]
            x = logZ(L, Jx, Jy, h, big(β))
            @test isapprox(Float64(mid(x)), ed; rtol = 1e-12)
        end
    end
end

@testset "closed form vs Arb referee" begin
    RR = Nemo.ArbField(1024)
    for (L, β) in ((16, 3 // 1), (64, 21 // 2), (128, 40 // 1)), (Jx, Jy, h) in PARAMS[1:3]
        x, ref = logZ(L, Jx, Jy, h, big(β)), arb_logZ(L, Jx, Jy, h, β)
        lo, hi = (RR(numerator(exact(e))) / RR(denominator(exact(e))) for e in (inf(x), sup(x)))
        @test lo <= ref && ref <= hi
        @test diam(x) < big(10)^-100
    end
end

@testset "domain and RP record" begin
    @test_throws ArgumentError logZ(8, -1 // 1, 1 // 2, 1 // 2, big(1) // 1)   # ordered phase
    @test_throws ArgumentError logZ(8, -1 // 1, -1 // 2, 2 // 1, big(1) // 1)  # J_x J_y > 0
    @test_throws ErrorException rp_record(1 // 2, -1 // 1, 3 // 2)  # needs the X ↔ Y rotation
    @test rp_record(-1 // 1, 1 // 4, 1 // 1)["weights"] == ["1//1", "1//4"]
end

@testset "stored certificates" begin
    for p in POINTS
        data = read_result(certificate_path(p), Dict{String,Any})
        b, c = check_point(data, p)
        @test 0 < c < b < exact_gap(p)
    end
end

@testset "tampering is rejected" begin
    p = POINTS[2]
    data = read_result(certificate_path(p), Dict{String,Any})
    bad = deepcopy(data)
    bad["seed"]["eps_t"] = string(parse_rational(data["seed"]["eps_t"]) / 2)
    @test_throws ErrorException check_point(bad, p)
    # With the bounds also forged to match, only the full replay of log Z catches it.
    ib, chain = bounds(parse_rational(bad["seed"]["eps_s"]),
        parse_rational(bad["seed"]["eps_t"]), p.n0, p.beta0)
    bad["seed"]["bound_exact_binary"] = string(RationalBig(ib.bound))
    bad["chain"] = chain
    @test check_point(bad, p) isa Tuple
    @test_throws ErrorException check_point(bad, p; full = true)
    @test_throws ErrorException check_point(data, POINTS[1])
end
