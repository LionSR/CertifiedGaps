"""
    IsingOpenSeed

Study `ising_open`: kernel-checked seed arithmetic for the longitudinal Ising chain with free
ends at g = 1/4, h = 1/20, the open-chain counterpart of `IsingSeed`. UNTRUSTED generator:
Float64 eigensolves only propose integer trial vectors and integer comparison weights; every
claim is re-evaluated by Lean's kernel in `lean/SeedKernel/IsingOpen/Generated*.lean`.

For each chain length n the generator writes

* trial vectors (one for n = 8, two for n = 2, 12) with their Gram and residual data, as for the
  rings of `IsingSeed`;
* for n = 2, 12, a cutoff `cut` with exactly two basis states below it, and positive integer
  weights `w` on the states above it with `5 ∑ w(flip) ≤ w(s) (20 E(s) - rho)`, the
  Collatz–Wielandt form of `H ≥ rho / 20` on the complement of those two states.
"""
module IsingOpenSeed

using LinearAlgebra, JSON, SHA
using ..Common: MODEL, LEAN, results_dir, export_all
using ..IsingSeed: SCALE, WIDTH, BIAS, pack, exactstr, exactint
BLAS.set_num_threads(1)
const OUT = results_dir("ising_open")
const LEANDIR = joinpath(LEAN, "SeedKernel", "IsingOpen")
"The generator, whose SHA-256 the certificate records (path relative to MODEL)."
const GENERATOR = "src/ising_open.jl"
"Chain lengths, the number of trial vectors, and whether comparison weights are needed."
const POINTS = [(2, 2, true), (8, 1, false), (12, 2, true)]
"Blocks of the 12-site checks: 64 blocks of 64 basis states."
const BLOCK = 64

"Integer diagonal of 20 H on the open chain: -20 ∑_{j<n-1} z_j z_{j+1} - ∑_j z_j."
function diagonal(n, s)
    z(j) = 1 - 2 * ((s >> j) & 1)
    -20 * sum((z(j) * z(j + 1) for j in 0:n-2); init=0) - sum(z(j) for j in 0:n-1)
end
"Apply the exact integer 20 H, without BLAS or a floating-point model."
apply_h(n, v) = [diagonal(n, s) * v[s+1] - 5sum(v[xor(s, 1 << j)+1] for j in 0:n-1)
                 for s in 0:(1<<n)-1]
"The weighted neighbour sum over flips that stay above the cutoff (Lean `flipSum`)."
flip_sum(n, w, cut, s) = sum((diagonal(n, xor(s, 1 << j)) < cut ? big(0) : w[xor(s, 1 << j)+1]
                              for j in 0:n-1); init=big(0))
"The comparison condition at state s (Lean `weightOK`)."
weight_ok(n, w, cut, rho, s) = diagonal(n, s) < cut ||
    (w[s+1] > 0 && 5flip_sum(n, w, cut, s) <= w[s+1] * (diagonal(n, s) - rho))

function float_h(n)
    d = 1 << n
    H = zeros(Float64, d, d)
    for s in 0:d-1
        H[s+1, s+1] = diagonal(n, s)
        for j in 0:n-1; H[xor(s, 1 << j)+1, s+1] = -5; end
    end
    H
end

"""
Positive integer weights on the states above the cutoff: the Perron vector of the compressed
`20 H` (rounded, scaled so that the largest is about 2^52), and the largest `rho` the integer
condition then allows, minus a margin of 1 (that is, 1/20 in units of H).
"""
function weights(n, cut)
    d = 1 << n
    keep = [s for s in 0:d-1 if diagonal(n, s) >= cut]
    Hc = float_h(n)[keep.+1, keep.+1]
    F = eigen(Symmetric(Hc), 1:1)
    p = abs.(F.vectors[:, 1])
    w = zeros(BigInt, d)
    for (i, s) in enumerate(keep)
        w[s+1] = max(big(1), round(BigInt, p[i] / maximum(p) * Float64(SCALE)))
    end
    # rho ≤ diagonal(s) - 5 flipSum / w(s) at every kept state
    bound = minimum(diagonal(n, s) - cld(5flip_sum(n, w, cut, s), w[s+1]) for s in keep)
    rho = bound - 1
    @assert all(s -> weight_ok(n, w, cut, rho, s), 0:d-1)
    w, rho
end

function proposal(n, nvec, weighted)
    d = 1 << n
    H = float_h(n)
    t = @elapsed F = eigen(Symmetric(H), 1:nvec)
    V = [round.(BigInt, F.vectors[:, j] .* Float64(SCALE)) for j in 1:nvec]
    lam = round.(BigInt, F.values .* Float64(SCALE))
    R = [SCALE .* apply_h(n, V[j]) .- lam[j] .* V[j] for j in 1:nvec]
    gram = [dot(V[i], V[j]) for i in 1:nvec, j in 1:nvec]
    res = [dot(r, r) for r in R]
    # All tests below are exact integers, even though their candidates came from LAPACK.
    @assert maximum(abs.(gram - SCALE^2 * I)) * big(10)^12 <= SCALE^2
    @assert sum(res) * big(10)^18 <= SCALE^4
    p = Dict{String,Any}("n" => n, "packed" => exactstr.(pack.(V)), "lambda" => exactstr.(lam),
        "gram" => exactstr.(vec(gram)), "residual_sq" => exactstr.(res),
        "generator_seconds" => t, "numerical_E" => F.values ./ 20)
    if weighted
        # The cutoff sits just above the second classical level: all spins up, all spins down.
        levels = sort(unique(diagonal(n, s) for s in 0:d-1))
        cut = levels[3]
        @assert count(s -> diagonal(n, s) < cut, 0:d-1) == 2
        w, rho = weights(n, cut)
        @assert lam[2] + SCALE < rho * SCALE  # the second eigenvalue lies below rho / 20
        p["cut_20H"] = exactstr(cut)
        p["rho_20H"] = exactstr(rho)
        p["weights"] = exactstr(pack(w))
    end
    p
end

lean_vec(n, i) = "u$(n)_$(i)"

"Lean source of `Generated.lean`: the proposals as constants and the small direct checks."
function generated_lean(cert)
    io = IOBuffer()
    println(io, "import SeedKernel.IsingOpen.Checker\n\nnamespace SeedKernel.IsingOpen\n")
    println(io, "open SeedKernel.Ising\n")
    println(io, "set_option maxRecDepth 1000000\nset_option maxHeartbeats 0\n")
    for p in cert["points"]
        n = p["n"]
        for i in eachindex(p["packed"])
            println(io, "def $(lean_vec(n, i)) : Nat := $(exactint(p["packed"][i]))")
            println(io, "def lam$(n)_$(i) : Int := $(exactint(p["lambda"][i]))")
        end
        if haskey(p, "weights")
            println(io, "def w$n : Nat := $(exactint(p["weights"]))")
            println(io, "def cut$n : Int := $(exactint(p["cut_20H"]))")
            println(io, "def rho$n : Int := $(exactint(p["rho_20H"]))")
            println(io, "theorem count$n : classicalCount $n cut$n = 2 := by decide +kernel")
        end
        n == 12 && continue  # checked block by block in GeneratedChunks.lean
        for i in eachindex(p["packed"])
            println(io, "theorem dot$(n)_$(i) : (dot $n $(lean_vec(n, i)) $(lean_vec(n, i)) - ",
                "scale ^ 2).natAbs * 10 ^ 12 ≤ (scale ^ 2).toNat := by decide +kernel")
        end
        res = join(["residualSq $n $(lean_vec(n, i)) lam$(n)_$(i)" for i in eachindex(p["packed"])],
            " + ")
        println(io, "theorem res$n : ($res) * 10 ^ 18 ≤ (scale ^ 4).toNat := by decide +kernel")
        haskey(p, "weights") && println(io,
            "theorem weights$n : weightBlock $n w$n cut$n rho$n 0 $(1 << n) = true := by decide +kernel")
    end
    println(io, "\nend SeedKernel.IsingOpen")
    String(take!(io))
end

"""
Lean source of `GeneratedChunks.lean`: the 12-site statistics and comparison conditions in 64
separately checked blocks of 64 states, which bound kernel memory, with their lists tied to a
canonical enumeration by `rfl`.
"""
function chunks_lean(cert)
    p = only(filter(q -> q["n"] == 12, cert["points"]))
    n = 12; d = 1 << n; nb = d ÷ BLOCK
    V = map(p["packed"]) do x
        u = exactint(x); [((u >> (WIDTH * s)) & ((big(1) << WIDTH) - 1)) - BIAS for s in 0:d-1]
    end
    lam = exactint.(p["lambda"])
    R = [SCALE .* apply_h(n, V[i]) .- lam[i] .* V[i] for i in 1:2]
    io = IOBuffer()
    println(io, "import SeedKernel.IsingOpen.Generated\n\nnamespace SeedKernel.IsingOpen\n")
    println(io, "open SeedKernel.Ising\n")
    println(io, "set_option maxRecDepth 1000000\nset_option maxHeartbeats 0\n")
    for i in 1:2, b in 0:nb-1
        ix = BLOCK*b+1:BLOCK*b+BLOCK
        println(io, "theorem stats12_$(i)_$b : blockStats 12 u12_$i lam12_$i $(BLOCK*b) $BLOCK = ",
            "($(sum(V[i][ix].^2)), $(sum(R[i][ix].^2))) := by decide +kernel")
    end
    for i in 1:2
        println(io, "def statBlocks12_$i : List (Int × Nat) := [")
        println(io, join(["  blockStats 12 u12_$i lam12_$i $(BLOCK*b) $BLOCK" for b in 0:nb-1],
            ",\n"), "]")
        println(io, "theorem statBlocks12_$(i)_eq : statBlocks12_$i = (List.range $nb).map ",
            "(fun b => blockStats 12 u12_$i lam12_$i ($BLOCK * b) $BLOCK) := by rfl")
        println(io, "theorem statSum12_$i : sumBlocks statBlocks12_$i = ",
            "($(sum(V[i].^2)), $(sum(R[i].^2))) := by")
        println(io, "  unfold statBlocks12_$i\n  rw [",
            join(["stats12_$(i)_$b" for b in 0:nb-1], ", "), "]\n  decide +kernel")
    end
    for b in 0:nb-1
        println(io, "theorem weights12_$b : weightBlock 12 w12 cut12 rho12 $(BLOCK*b) $BLOCK = true := ",
            "by decide +kernel")
    end
    println(io, "def weightBlocks12 : List Bool := [")
    println(io, join(["  weightBlock 12 w12 cut12 rho12 $(BLOCK*b) $BLOCK" for b in 0:nb-1], ",\n"), "]")
    println(io, "theorem weightBlocks12_eq : weightBlocks12 = (List.range $nb).map ",
        "(fun b => weightBlock 12 w12 cut12 rho12 ($BLOCK * b) $BLOCK) := by rfl")
    println(io, "theorem weightBlocks12_all : weightBlocks12.all id = true := by")
    println(io, "  unfold weightBlocks12\n  rw [", join(["weights12_$b" for b in 0:nb-1], ", "),
        "]\n  decide +kernel")
    println(io, "\nend SeedKernel.IsingOpen")
    String(take!(io))
end

"""
    certify_ising_open()

Propose vectors and weights for the open chains n = 2, 8, 12 (Float64 LAPACK, untrusted), round
them to exact integers, and write `results/ising_open/certificate.json` and
`lean/SeedKernel/IsingOpen/Generated{,Chunks}.lean`. About 20 s.
"""
function certify_ising_open()
    cert = Dict("scale" => exactstr(SCALE), "width" => WIDTH, "bias" => exactstr(BIAS),
        "points" => [proposal(p...) for p in POINTS],
        "source_hashes" => Dict(GENERATOR => bytes2hex(sha256(read(joinpath(MODEL, GENERATOR))))),
        "run_metadata" => Dict("julia" => string(VERSION), "blas_threads" => 1,
            "JSON" => string(Base.pkgversion(JSON)), "blas" => string(BLAS.get_config()),
            "generator" => "untrusted Float64 LAPACK plus exact BigInt arithmetic"))
    mkpath(OUT)
    open(joinpath(OUT, "certificate.json"), "w") do io; JSON.print(io, cert); end
    # Round-trip through the file so that the Lean files are generated from what is stored.
    stored = JSON.parsefile(joinpath(OUT, "certificate.json"))
    @assert verify_certificate(stored)
    mkpath(LEANDIR)
    write(joinpath(LEANDIR, "Generated.lean"), generated_lean(stored))
    write(joinpath(LEANDIR, "GeneratedChunks.lean"), chunks_lean(stored))
    println("Generated exact integer proposals and weights for the open chains n = 2, 8, 12.")
end

"Unpack a stored vector of `2^n` biased 56-bit digits."
unpack(x, n) = (u = exactint(x); [((u >> (WIDTH * s)) & ((big(1) << WIDTH) - 1)) - BIAS
                                  for s in 0:(1<<n)-1])

"Exact Julia replay of every integer claim of the certificate (a diagnostic)."
function verify_certificate(c)
    @assert exactint(c["scale"]) == SCALE
    @assert c["width"] == WIDTH && exactint(c["bias"]) == BIAS
    @assert [(p["n"], length(p["packed"]), haskey(p, "weights")) for p in c["points"]] == POINTS
    for p in c["points"]
        n = p["n"]; d = 1 << n; k = length(p["packed"])
        V = [unpack(x, n) for x in p["packed"]]
        @assert exactstr.(pack.(V)) == p["packed"] # rejects unconsumed high digits
        lam = exactint.(p["lambda"])
        G = [sum(V[i] .* V[j]) for i in 1:k, j in 1:k]
        @assert exactstr.(vec(G)) == p["gram"]
        R = [SCALE .* apply_h(n, V[i]) .- lam[i] .* V[i] for i in 1:k]
        rs = [sum(r .^ 2) for r in R]
        @assert exactstr.(rs) == p["residual_sq"]
        @assert maximum(abs.(G - SCALE^2 * I)) * big(10)^12 <= SCALE^2
        @assert sum(rs) * big(10)^18 <= SCALE^4
        if haskey(p, "weights")
            cut = exactint(p["cut_20H"]); rho = exactint(p["rho_20H"])
            w = unpack(p["weights"], n)
            @assert exactstr(pack(w)) == p["weights"]
            @assert count(s -> diagonal(n, s) < cut, 0:d-1) == 2
            @assert all(s -> weight_ok(n, w, cut, rho, s), 0:d-1)
            @assert lam[2] + SCALE < rho * SCALE
        end
    end
    true
end

"""
    check_ising_open()

Exact replay of `results/ising_open/certificate.json`, the recorded generator hash, and byte
equality of `Generated.lean`/`GeneratedChunks.lean` with the text generated from the stored
certificate. The Lean kernel itself is run by `Common.build_seed_kernel`. Throws on mismatch.
"""
function check_ising_open()
    c = JSON.parsefile(joinpath(OUT, "certificate.json"))
    @assert verify_certificate(c)
    @assert c["source_hashes"] ==
        Dict(GENERATOR => bytes2hex(sha256(read(joinpath(MODEL, GENERATOR)))))
    generated_lean(c) == read(joinpath(LEANDIR, "Generated.lean"), String) ||
        error("lean/SeedKernel/IsingOpen/Generated.lean differs from the certificate")
    chunks_lean(c) == read(joinpath(LEANDIR, "GeneratedChunks.lean"), String) ||
        error("lean/SeedKernel/IsingOpen/GeneratedChunks.lean differs from the certificate")
    println("ising_open: exact Julia replay, generator hash and generated Lean files match")
    true
end

export_all(@__MODULE__)

end
