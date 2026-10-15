using Test
include("../src/TransferCertificate.jl")
const TC = TransferCertificate

@testset "Transfer inequalities and Weyl transport" begin
    # Exact probability vectors test the general inequality, including equality at rank two.
    for weights in ([1, 0], [9, 1], [6, 3, 1], [1, 1, 1], [17, 2, 1, 1])
        p = big.(weights) .// sum(weights)
        purity = sum(p.^2)
        defect = 1 - purity
        squared = 1 - sum(p.^4) / purity^2
        @test squared <= defect^2 / (2 * purity^2)
        @test squared <= TC.Q(TC.sup(TC.square_defect(TC.box(defect))))
    end
    @test TC.sup(TC.square_defect(TC.box(1))) == 1
    @test_throws ArgumentError TC.square_defect(TC.box(-1//10))

    # The tail includes all omitted multiplicities even when their individual radii differ.
    spec = (; E = [0.0, 1.0, 2.0, 3.0], R = [0.01, 0.02, 0.04, 0.03], L = 2)
    full = TC.logz(spec, 3//2; keep = 4)
    truncated = TC.logz(spec, 3//2; keep = 2)
    @test TC.inf(truncated) <= TC.inf(full) <= TC.sup(full) <= TC.sup(truncated)
    @test_throws ArgumentError TC.logz(spec, 0//1)
    @test_throws ArgumentError TC.logz(spec, 1//1; keep = 0)

    # At L=2 this compares independently enclosed spectra at different rational D.
    quiet(D) = redirect_stdout(() -> TC.compute_spectrum(2, D; rigorous = true), devnull)
    specs = [quiet(D) for D in (3//1, 31//10)]
    z = [TC.logz((; E = s["energies"], R = s["radii"]), 2//1) for s in specs]
    moved = TC.transported_logz(z[1], 2, 2//1, 3//1, 31//10, 31//10)
    shifted = z[2] + TC.box(2//10)
    @test TC.inf(moved) < TC.inf(shifted) < TC.sup(shifted) < TC.sup(moved)
    @test_throws ArgumentError TC.transported_logz(z[1], 2, 2//1, 3//1, 4//1, 3//1)

    # Coupling is an exact cancellation of partition functions, independent of a model.
    a, b, c, d = big.([11, 17, 23, 31]) .// 1
    ps, ps2, pt, pt2 = b/a^2, d/c^2, c/a^2, d/b^2
    @test pt2 == ps2 * pt^2 / ps^2

    # A level from a later block can cross the first two centers; use a common radius.
    crossing = (; E = [0.0, 1.0, 1.1], R = [0.0, 0.0, 1.0], L = 2, D = 4//1)
    @test TC.inf(TC.finite_gap(crossing, 4//1, 4//1)) <= 1//10
end

@testset "Certificate replay rejects altered hypotheses" begin
    a, b = TC.reference_spectrum(4//1, 6), TC.reference_spectrum(4//1, 12)
    row = TC.evaluate(a, b, 29//10; max_lifts = 1)[end]
    @test row.n0 == 12 && row.result.bound > 8//10
    data = TC.record(row, b, "mutation-test")
    @test TC.check_record(data, a, b).result.bound == row.result.bound
    for (key, value) in (("n0", "7"), ("bound", "2//1"), ("eps_s", "0//1"),
        ("reference_D", "3//1"), ("D_interval", ["4//1", "3//1"]), ("logZ", Any[]))
        changed = deepcopy(data)
        changed[key] = value
        @test_throws ErrorException TC.check_record(changed, a, b)
    end
end
