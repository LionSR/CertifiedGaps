# Fast tests of the cluster-chain certificate.
#
#     julia --project=code/julia code/models/cluster_chain/test/runtests.jl
using Test
include(joinpath(@__DIR__, "..", "src", "common.jl"))
import Nemo

# Dense H on N sites in the σᶻ product basis (bit j of the index is site j: 0 ↔ Z = +1).
function dense_cluster(N, g)
    H = zeros(2^N, 2^N)
    for b in 0:(2^N-1)
        z(j) = iseven(b >> mod(j, N)) ? 1 : -1
        for j in 0:(N-1)
            H[(b⊻(1<<j))+1, b+1] -= z(j - 1) * z(j + 1) + g
        end
    end
    return H
end

# Dense ring Hamiltonian of a 4-dimensional-cell RingModel on L cells, in the same site basis:
# cell c holds sites (2c, 2c+1) and its local index is 2·bit(2c) + bit(2c+1), the first
# Kronecker factor being site 2c.
function dense_from_cells(m, L)
    N = 2L
    cell(b, c) = 2 * ((b >> 2c) & 1) + ((b >> (2c + 1)) & 1)
    setcell(b, c, s) = (b & ~(3 << 2c)) | ((s >> 1) << 2c) | ((s & 1) << (2c + 1))
    H = zeros(ComplexF64, 2^N, 2^N)
    for b in 0:(2^N-1), c in 0:(L-1)
        c2 = mod(c + 1, L)
        s1, s2 = cell(b, c), cell(b, c2)
        for t1 in 0:3, t2 in 0:3
            v = m.bond[4t1+t2+1, 4s1+s2+1]
            iszero(v) || (H[setcell(setcell(b, c, t1), c2, t2)+1, b+1] += ComplexF64(v))
        end
        for t in 0:3
            v = m.onsite[t+1, s1+1]
            iszero(v) || (H[setcell(b, c, t)+1, b+1] += ComplexF64(v))
        end
    end
    return H
end

# Independent evaluation of the closed form in Arb balls (Nemo), straight from the four
# products. Rationals are built from their integers: Nemo converts a Julia Rational through a
# BigFloat at the global precision.
function arb_logZ(L, g, β; bits = 1024)
    RR = Nemo.ArbField(bits)
    q(x) = RR(numerator(x)) / RR(denominator(x))
    G, b, N = q(g), q(β), 2L
    e(m, s) = 2 * sqrt(1 + G^2 - 2G * cos(4 * Nemo.const_pi(RR) * (m + q(s)) / N))
    prods(s) = (prod(2cosh(b * e(m, s) / 2) for m in 0:(N-1)),
        prod(2sinh(b * e(m, s) / 2) for m in 0:(N-1)))
    (cn, sn), (cr, sr) = prods(1 // 2), prods(0 // 1)
    return log((cn + sn + cr - sr) / 2)
end

@testset "cell model reproduces the site Hamiltonian" begin
    for L in 2:4, g in (1 // 2, 2 // 1)
        @test dense_from_cells(physical_model(g), L) ≈ dense_cluster(2L, Float64(g))
        # The rotated model has the same spectrum (U is an exact frame change).
        @test eigvals(Hermitian(dense_from_cells(rotated_model(g), L))) ≈
              eigvals(Symmetric(dense_cluster(2L, Float64(g))))
    end
end

@testset "closed form vs exact diagonalization" begin
    for L in 2:5, g in (1 // 4, 1 // 2, 2 // 1, 3 // 1)
        ev = eigvals(Symmetric(dense_cluster(2L, Float64(g))))
        # Unique ground state with gap ≥ 2|1 − g| on these rings.
        @test ev[2] - ev[1] > 2abs(1 - g) - 1e-9
        for β in (1 // 3, 2 // 1, 7 // 1)
            ed = log(sum(exp.(-Float64(β) .* (ev .- ev[1])))) - Float64(β) * ev[1]
            x = logZ(L, g, big(β))
            @test isapprox(Float64(mid(x)), ed; rtol = 1e-12)
        end
    end
end

@testset "closed form vs Arb referee" begin
    RR = Nemo.ArbField(1024)
    for (L, g, β) in ((8, 1 // 2, 3 // 1), (32, 1 // 2, 21 // 2), (64, 1 // 4, 40 // 1))
        x, ref = logZ(L, g, big(β)), arb_logZ(L, g, β)
        lo, hi = (RR(numerator(exact(e))) / RR(denominator(exact(e))) for e in (inf(x), sup(x)))
        @test lo <= ref && ref <= hi
        @test diam(x) < big(10)^-100
    end
end

@testset "domain and RP record" begin
    @test_throws ArgumentError logZ(4, 1 // 1, big(1) // 1)
    @test_throws ArgumentError mode_energies(5, 1 // 2, :ns)
    r = rp_record(1 // 2)
    @test !isempty(r["factors"]) && all(w -> parse_rational(w) > 0, r["weights"])
    # The physical product basis is not reflection positive.
    @test_throws ArgumentError rp_certificate(physical_model(1 // 2))
end

@testset "stored certificate and tampering" begin
    p = POINTS[1]
    data = read_result(certificate_path(p.g), Dict{String,Any})
    b, c = check_point(data, p.g, p.n0, p.beta0)
    @test 0 < c < b < 2abs(1 - p.g)
    bad = deepcopy(data)
    bad["seed"]["eps_s"] = string(parse_rational(data["seed"]["eps_s"]) / 2)
    @test_throws ErrorException check_point(bad, p.g, p.n0, p.beta0)
    ib, chain = bounds(parse_rational(bad["seed"]["eps_s"]),
        parse_rational(bad["seed"]["eps_t"]), p.n0, p.beta0)
    bad["seed"]["bound_exact_binary"] = string(RationalBig(ib.bound))
    bad["chain"] = chain
    @test check_point(bad, p.g, p.n0, p.beta0) isa Tuple
    @test_throws ErrorException check_point(bad, p.g, p.n0, p.beta0; full = true)
    bad2 = deepcopy(data)
    bad2["reflection_positivity"]["weights"][1] = "2"
    @test_throws ErrorException check_point(bad2, p.g, p.n0, p.beta0)
end
