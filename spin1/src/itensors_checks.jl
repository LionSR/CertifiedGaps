"Numerical ITensors controls for the spin-1 certificates; never a certificate."
module ITensorsChecks
using ITensors, ITensorMPS, LinearAlgebra, Random, Dates
using ThermalMPOGap.Purification: purified_log_partition
import ..Spin1Certificates as SC

function bonds(n)
    n >= 4 || throw(ArgumentError("periodic spin-1 chains require n >= 4"))
    return [(j, mod1(j + 1, n)) for j in 1:n]
end

function hamiltonian(sites; twist = 0)
    twist in (0, 1, 2) || throw(ArgumentError("twist must be 0, 1 or 2"))
    os = OpSum()
    for (j, k) in bonds(length(sites))
        phase = k != 1 || twist == 0 ? 1.0 : twist == 1 ? -1.0 : cis(2pi / 3)
        os += 1.0, "Sz", j, "Sz", k
        os += phase / 2, "S+", j, "S-", k
        os += conj(phase) / 2, "S-", j, "S+", k
    end
    return MPO(os, sites)
end

function bond_operator(sites, j, k)
    return op("Sz", sites, j) * op("Sz", sites, k) +
           (op("S+", sites, j) * op("S-", sites, k) + op("S-", sites, j) * op("S+", sites, k)) / 2
end

"Numerical second-order thermal MPO control, not a certificate (core purification)."
function thermal_logz(n, beta; dbeta = 0.025, maxdim = 81, cutoff = 1e-14)
    isfinite(beta) && beta >= 0 || throw(ArgumentError("beta must be finite and nonnegative"))
    isfinite(dbeta) && dbeta > 0 || throw(ArgumentError("dbeta must be finite and positive"))
    maxdim > 0 && isfinite(cutoff) && cutoff >= 0 ||
        throw(ArgumentError("invalid truncation parameters"))
    steps = round(Int, beta / dbeta)
    isapprox(steps * dbeta, beta; atol = 1e-12, rtol = 1e-12) ||
        throw(ArgumentError("beta must be a multiple of dbeta"))
    sites = siteinds("S=1", n; conserve_qns = true)
    halfgate = (sites, i, j, tau) -> exp(-tau * bond_operator(sites, i, j))
    result = purified_log_partition(sites, bonds(n), halfgate, [beta]; dbeta, maxdim, cutoff)
    return result[Float64(beta)]["logZ"]
end

function dense_physical_reference(n, twist)
    h = SC.physical_hamiltonian(Val(n), twist)
    result = zeros(ComplexF64, 3^n, 3^n)
    index(word) = 1 + sum((1 - word[j]) * 3^(j - 1) for j in 1:n)
    omega = cis(2pi / 3)
    for ((x, y), (a, b)) in h
        result[index(x), index(y)] = a + b * omega
    end
    return result
end

function run_itensors_checks(output)
    BLAS.set_num_threads(1)
    started = time()
    sites = siteinds("S=1", 4)
    matrix_errors = Float64[]
    traces = Any[]
    spectrum = Float64[]
    for twist in 0:2
        H = hamiltonian(sites; twist)
        dense = reshape(Array(prod(H), prime.(sites)..., dag.(sites)...), 81, 81)
        physical = dense_physical_reference(4, twist)
        error = maximum(abs, dense - physical)
        @assert error < 1e-12 && isapprox(dense, dense'; atol = 1e-12)
        push!(matrix_errors, error)
        eigenvalues = eigvals(Hermitian(dense))
        twist == 0 && (spectrum = eigenvalues)
        if twist <= 1
            beta, a = 49 / 4, 700741 / 500000
            trace = sum(exp.(-beta .* (eigenvalues .+ 4a)))
            expected = only(
                filter(
                    x -> SC.integer(x["n"]) == 4 && SC.integer(x["tw"]) == twist,
                    SC.read_exact_json(
                        joinpath(SC.SOURCE, "evidence/source1/thermal-run/certificates.json"),
                    ),
                ),
            )
            lower = Float64(SC.parse_rational(expected["lower"]))
            upper = Float64(SC.parse_rational(expected["upper"]))
            @assert lower < trace < upper
            push!(
                traces,
                Dict("twist" => twist, "trace" => trace, "within_rigorous_interval" => true),
            )
        end
    end
    @assert isapprox(spectrum[1], -6; atol = 1e-12)
    @assert isapprox(spectrum[2] - spectrum[1], 1; atol = 1e-12)
    _, A = SC.trial_matrices()
    amplitudes = zeros(Float64, 3, 3, 3, 3)
    for word in Iterators.product(ntuple(_ -> -1:1, 4)...)
        amplitudes[map(s->2-s, word)...] = Float64(tr(prod(A[s] for s in word)))
    end
    amplitudes /= norm(amplitudes)
    psi = MPS(ITensor(amplitudes, sites...), sites; cutoff = 0.0)
    trial_energy = real(inner(psi', hamiltonian(sites), psi)) / real(inner(psi, psi))
    exact_energy = Float64(
        SC.frac(
            parse(BigInt, "-137851234864097615402880"),
            parse(BigInt, "22975205960525898253776"),
        ),
    )
    @assert isapprox(trial_energy, exact_energy; atol = 1e-11)
    qsites = siteinds("S=1", 4; conserve_qns = true)
    Random.seed!(1729)
    initial = random_mps(qsites, ["Up", "Dn", "Up", "Dn"]; linkdims = 9)
    dmrg_energy, _ = dmrg(
        hamiltonian(qsites),
        initial;
        nsweeps = 6,
        maxdim = [16, 27],
        cutoff = 1e-13,
        outputlevel = 0,
    )
    @assert isapprox(dmrg_energy, -6; atol = 1e-9)
    beta = 0.1
    exact_logz = log(sum(exp.(-beta .* spectrum)))
    coarse = thermal_logz(4, beta; dbeta = 0.05)
    fine = thermal_logz(4, beta; dbeta = 0.025)
    @assert abs(fine - exact_logz) < abs(coarse - exact_logz) / 3
    @assert abs(fine - exact_logz) < 1e-5
    @assert isapprox(thermal_logz(4, 0.0), 4log(3); atol = 1e-12)
    println(
        "PASS ITensors: three twisted MPOs, trial MPS energy, DMRG energy -6, thermal MPO convergence",
    )
    return SC.write_result(
        output,
        Dict(
            "status" => "numerical-crosscheck-passed",
            "certified" => false,
            "julia" => string(VERSION),
            "ITensors" => string(pkgversion(ITensors)),
            "ITensorMPS" => string(pkgversion(ITensorMPS)),
            "matrix_max_errors" => matrix_errors,
            "thermal_traces" => traces,
            "trial_energy" => trial_energy,
            "exact_trial_energy" => exact_energy,
            "dmrg_energy" => dmrg_energy,
            "thermal_mpo" => Dict(
                "n" => 4,
                "beta" => beta,
                "exact_logZ" => exact_logz,
                "coarse_dbeta" => 0.05,
                "coarse_logZ" => coarse,
                "fine_dbeta" => 0.025,
                "fine_logZ" => fine,
                "fine_absolute_error" => abs(fine - exact_logz),
            ),
            "wall_seconds" => time() - started,
            "completed_utc" => string(now(UTC)),
        ),
    )
end
end
