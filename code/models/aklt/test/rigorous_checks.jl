using Test, LinearAlgebra

const AKLTExact = Rational{BigInt}

"""Maximal-spin projector in the rational binomial similarity gauge."""
function aklt_rational_projector()
    weights = BigInt[binomial(3, a) for a in 0:3]
    p = zeros(AKLTExact, 16, 16)
    for a in 0:3, b in 0:3, c in 0:3, d in 0:3
        if a + b == c + d
            p[4a + b + 1, 4c + d + 1] = weights[c + 1] * weights[d + 1] //
                binomial(big(6), a + b)
        end
    end
    return p, weights
end

"""Matrix of A -> Tr_2[p (I tensor A)], with column-major operator coordinates."""
function aklt_frame_matrix(p)
    k = zeros(AKLTExact, 16, 16)
    for i in 0:3, j in 0:3, l in 0:3, m in 0:3
        k[i + 4m + 1, l + 4j + 1] = p[4i + j + 1, 4m + l + 1]
    end
    return k
end

@testset "exact AKLT coherent-frame and RP obstruction" begin
    p, weights = aklt_rational_projector()
    i4 = Matrix{AKLTExact}(I, 4, 4)
    i16 = Matrix{AKLTExact}(I, 16, 16)
    w = Diagonal(AKLTExact.(weights))
    @test p * p == p
    @test tr(p) == 7
    @test kron(w, w) * p == transpose(p) * kron(w, w)

    # Independent construction from the total-spin polynomial, in the same gauge.
    sz = Diagonal(AKLTExact[3//2, 1//2, -1//2, -3//2])
    sp = diagm(1 => AKLTExact[3, 2, 1])
    sm = diagm(-1 => AKLTExact[1, 2, 3])
    x = kron(sz, sz) + (kron(sp, sm) + kron(sm, sp)) / 2
    polynomial = (x + 15i16 / 4) * (x + 11i16 / 4) * (x + 3i16 / 4) / 90
    @test p == polynomial

    k = aklt_frame_matrix(p)
    # The gauge is not unitary. This metric is the physical Hilbert-Schmidt metric.
    metric = Diagonal(AKLTExact[weights[a] // weights[b] for b in 1:4 for a in 1:4])
    @test metric * k == transpose(k) * metric
    @test k * vec(i4) == (7//4) * vec(i4)
    eigenvalues = AKLTExact[7//4, 21//20, 7//20, 1//20]
    multiplicities = [1, 3, 5, 7]
    # Independent exact integration of the four Funk-Hecke formulas.
    overlap = AKLTExact[1, 3, 3, 1] / 8
    legendre = [AKLTExact[1], AKLTExact[0, 1],
        AKLTExact[-1//2, 0, 3//2], AKLTExact[0, -3//2, 0, 5//2]]
    for (poly, lambda) in zip(legendre, eigenvalues)
        integral = zero(AKLTExact)
        for a in eachindex(overlap), b in eachindex(poly)
            degree = a + b - 2
            if iseven(degree)
                integral += overlap[a] * poly[b] * (big(2) // (degree + 1))
            end
        end
        @test (7//2) * integral == lambda
    end
    projectors = Matrix{AKLTExact}[]
    for (lambda, multiplicity) in zip(eigenvalues, multiplicities)
        e = copy(i16)
        for mu in eigenvalues
            mu == lambda && continue
            e = e * (k - mu * i16) / (lambda - mu)
        end
        @test e * e == e
        @test k * e == lambda * e
        @test metric * e == transpose(e) * metric
        @test tr(e) == multiplicity
        push!(projectors, e)
    end
    @test sum(projectors) == i16
    @test sum(eigenvalues .* projectors) == k
    @test minimum(eigenvalues) == 1//20
    for a in 1:4, b in (a + 1):4
        @test iszero(projectors[a] * projectors[b])
    end

    q = Diagonal(AKLTExact[1, -1, -1, 1])
    @test tr(q) == 0
    @test tr(q * q) == 4
    @test tr(p * kron(q, q)) == 7//5
    @test k * vec(Matrix(q)) == (7//20) * vec(Matrix(q))

    # Three-site trace checks the cancellation of a term not containing both witnesses.
    witness = kron(q / 2, q / 2, i4)
    direct = kron(p, i4)
    other = kron(i4, p)
    @test tr(other * witness) == 0
    @test -tr((direct + other) * witness) / 64 == -7//320
    @test (big(1)//80) / 8 == 1//640
    @test 1//320 - 1//640 == 1//640

    # Deliberate corruption must fail the exact identities, rather than pass a tolerance.
    bad = copy(p)
    bad[1, 1] += 1//1000
    @test bad * bad != bad
    badk = copy(k)
    badk[1, 1] += 1//1000
    @test badk * vec(i4) != (7//4) * vec(i4)
end

@testset "exact trajectory-envelope constants" begin
    rmax, etamax = big(1)//10, big(1)//40
    coefficient = (4 + 6rmax^2 + 4rmax^4 + rmax^6) / (1 - etamax)
    @test coefficient == 1353467//325000
    @test coefficient < 5
    @test inv(1 - etamax) < 2
    @test 5rmax^2 + 2etamax == rmax

    # exp(eta) <= 1/(1-eta); all subsequent arithmetic is rational and exact.
    for q in AKLTExact[1//100, 1//10, 1//2, 9//10], j in 0:6
        r = q^(2^j) / 10
        eta = q^(2^(j + 1)) / 40
        majorant = (1 + r^2)^4 / (1 - eta) - 1
        @test majorant <= 5r^2 + 2eta
        @test 5r^2 + 2eta == q^(2^(j + 1)) / 10
    end

    # Finite-product identity behind the weighted-log criterion; no log evaluation.
    c, power, r0 = big(3)//1, 2, big(1)//20
    r = r0
    for j in 0:8
        expected = (c * 2^power * r0)^(2^j) / (c * big(2)^(power * (j + 1)))
        @test r == expected
        @test sum((big(1)//2)^(k + 1) for k in 0:(j - 1); init=0//1) ==
            1 - (big(1)//2)^j
        @test sum(k * (big(1)//2)^(k + 1) for k in 0:(j - 1); init=0//1) ==
            1 - (j + 1) * (big(1)//2)^j
        r = c * big(2)^(power * j) * r^2
    end
end
