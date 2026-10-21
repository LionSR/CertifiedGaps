# Stage-2 inequalities on exact data, Weyl transport on independently verified spectra, and
# rejection of altered certificates.
using Test, IntervalArithmetic
using CertifiedGaps.Results: read_exact_json, write_result
using CertifiedGaps.VerifiedSpectra: verified_ring_spectrum, BlockSpectrum
isdefined(@__MODULE__, :TransferCertificate) ||
    include(joinpath(dirname(@__DIR__), "src", "TransferCertificate.jl"))
const TC = TransferCertificate

@testset "Transfer inequalities and Weyl transport" begin
    # Exact probability vectors test the general inequality, including equality at rank two.
    for weights in ([1, 0], [9, 1], [6, 3, 1], [1, 1, 1], [17, 2, 1, 1])
        p = big.(weights) .// sum(weights)
        purity = sum(p .^ 2)
        defect = 1 - purity
        squared = 1 - sum(p .^ 4) / purity^2
        @test squared <= defect^2 / (2 * purity^2)
        @test squared <= TC.Q(sup(TC.square_defect(TC.box(defect))))
    end
    @test sup(TC.square_defect(TC.box(1))) == 1
    @test_throws ArgumentError TC.square_defect(TC.box(-1 // 10))

    # At L = 2 this compares independently verified spectra at different rational D.
    spectra = [(; D, L = 2, spectrum = verified_ring_spectrum(TC.ring_model(D), 2))
               for D in (3 // 1, 31 // 10)]
    z = [TC.logz(s, 2 // 1) for s in spectra]
    moved = TC.transported_logz(z[1], 2, 2 // 1, 3 // 1, 31 // 10, 31 // 10)
    shifted = z[2] + TC.box(2 // 10)                 # the centering shift L δD / 2 · β
    @test inf(moved) < inf(shifted) < sup(shifted) < sup(moved)
    @test_throws ArgumentError TC.transported_logz(z[1], 2, 2 // 1, 3 // 1, 4 // 1, 3 // 1)
    gap = TC.finite_gap(spectra[1], 3 // 1, 31 // 10)
    exact = TC.finite_gap(spectra[2], 31 // 10, 31 // 10)
    @test gap[1] <= exact[1] <= exact[2] <= gap[2]

    # Coupling is an exact cancellation of partition functions, independent of a model.
    a, b, c, d = big.([11, 17, 23, 31]) .// 1
    ps, ps2, pt, pt2 = b / a^2, d / c^2, c / a^2, d / b^2
    @test pt2 == ps2 * pt^2 / ps^2
end

@testset "Certificate replay rejects altered hypotheses" begin
    a, b = TC.reference_spectrum(4 // 1, 6), TC.reference_spectrum(4 // 1, 12)
    row = TC.evaluate(a, b, 29 // 10; max_lifts = 1)[end]
    @test row.n0 == 12 && row.result.bound > 8 // 10
    data = TC.record(row, b, "mutation-test")
    @test TC.check_record(data, a, b).result.bound == row.result.bound
    # Weakened values ("1//10", "99//100") must fail too: the replay requires equality.
    for (key, value) in (("n0", "7"), ("bound", "2//1"), ("bound", "1//10"), ("eps_s", "0//1"),
        ("eps_t", "99//100"), ("reference_D", "3//1"), ("D_interval", ["4//1", "3//1"]),
        ("logZ", Any[]), ("finite_gap_L12", ["3//1", "4//1"]))
        changed = deepcopy(data)
        changed[key] = value
        @test_throws ErrorException TC.check_record(changed, a, b)
    end
end

@testset "Certificate replay rejects a truncated case list or a moved interval" begin
    data = read_exact_json(joinpath(TC.ROOT, "results", "certificate.json"))
    truncated = deepcopy(data)
    filter!(case -> case["id"] != "transport-D1603//500", truncated["cases"])
    path = tempname() * ".json"
    write_result(path, truncated)
    @test_throws ErrorException TC.check_certificate(path)
    # A self-consistent case on a narrower interval than the advertised one is refused too.
    a, b = TC.reference_spectrum(1603 // 500, 6), TC.reference_spectrum(1603 // 500, 12)
    row = TC.evaluate(a, b, 381 // 125; lo = 32059999 // 10^7, hi = 1603 // 500,
        max_lifts = 0)[end]
    moved = deepcopy(data)
    index = findfirst(case -> case["id"] == "transport-D1603//500", moved["cases"])
    moved["cases"][index] = TC.record(row, b, "transport-D1603//500")
    @test TC.check_record(moved["cases"][index], a, b) !== nothing
    write_result(path, moved)
    @test_throws "advertised" TC.check_certificate(path)
end

@testset "Curve data comes from the committed certificate" begin
    curve = read_exact_json(joinpath(TC.ROOT, "results", "transfer_curve.json"))
    cases = read_exact_json(joinpath(TC.ROOT, "results", "certificate.json"))["cases"]
    @test length(curve["certified"]) == 2
    for point in curve["certified"]
        case = only(filter(c -> c["id"] == point["id"], cases))
        @test point["D_interval"] == case["D_interval"]
        @test point["bound"] == case["bound"]
        @test point["gap_L12_interval"] == case["finite_gap_L12"]
        @test point["ratio_interval"] == case["ratio_to_gap_L12"]
    end
end
