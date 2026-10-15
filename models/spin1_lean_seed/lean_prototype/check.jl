# Check the spin-one packed Lean prototype. Julia never decides the truth of a Lean statement;
# it regenerates the witness recipes, compares them with the stored inventory, and runs Lean.
#
#     julia --project=. check.jl          # hashes, regeneration, clean `lake build`, axioms
#     julia --project=. check.jl --full   # also the axiom audit and every src/Bench file
#
# Set LAKE_BIN if `lake` is not on PATH. Fails loudly (throws) on any mismatch.
using CertifiedGaps.Results: read_exact_json

include(joinpath(@__DIR__, "scripts", "generate.jl"))
include(joinpath(@__DIR__, "scripts", "witnesses.jl"))

const LAKE = get(ENV, "LAKE_BIN", "lake")
const SRC = joinpath(@__DIR__, "src")
const STANDARD_AXIOMS = Set(["propext", "Classical.choice", "Quot.sound"])

"""Run `cmd` in src/, return its combined output; error (with the output) if it fails."""
function lean_output(cmd)
    buffer = IOBuffer()
    process = run(pipeline(ignorestatus(Cmd(cmd; dir=SRC)); stdout=buffer, stderr=buffer))
    text = String(take!(buffer))
    process.exitcode == 0 || error("$(cmd) failed:\n$text")
    occursin(r"\berror\b|declaration uses 'sorry'", text) && error("$(cmd) reported:\n$text")
    return text
end

"""Every `#print axioms` line in `text` must name only propext, Classical.choice, Quot.sound."""
function check_axioms(text; expected)
    seen = String[]
    for m in eachmatch(r"'([^']+)' depends on axioms: \[([^\]]*)\]", text)
        axioms = Set(strip.(split(m[2], ',')))
        issubset(axioms, STANDARD_AXIOMS) || error("$(m[1]) uses nonstandard axioms $(axioms)")
        push!(seen, m[1])
    end
    missing_names = setdiff(expected, seen)
    isempty(missing_names) || error("no axiom report for $(missing_names)")
    return seen
end

function main(full)
    stored = read_exact_json(joinpath(@__DIR__, "results", "witnesses.json"))
    recorded = Dict(string(k) => string(v) for (k, v) in stored["run_metadata"]["source_sha256"])
    current = hashed_sources(@__DIR__)
    recorded == current || error("source hashes differ from results/witnesses.json:\n" *
        join(sort(collect(symdiff(Set(recorded), Set(current)))), "\n"))
    println("source hashes match ($(length(current)) files)")

    for (row, (n, q)) in zip(stored["witnesses"], WITNESS_SIZES)
        fresh, text = generate(n, q; write=false)
        _, _, width, num, den, epsnum, epsden, g, a, c = fresh
        expected = Dict("n" => "$n", "degree" => "$q", "lane_bits" => "$width",
            "ritz_upper" => "-$num//$den", "residual_upper" => "$epsnum//$epsden",
            "G" => "$g", "A_filter" => "$a", "C_filter" => "$c")
        Dict(string(k) => string(v) for (k, v) in row) == expected ||
            error("H$n witness differs from regeneration")
        text == read(joinpath(SRC, "SeedPrototype", "H$n.lean"), String) ||
            error("src/SeedPrototype/H$n.lean differs from the generator's recipe")
        println("H$n: Ritz <= -$num/$den, residual <= $epsnum/$epsden, $width-bit lanes")
    end

    # Clean rebuild of the default target so that every theorem is re-checked by the kernel.
    build = ""
    t = @elapsed begin
        lean_output(`$LAKE clean`)
        build = lean_output(`$LAKE build`)
    end
    names = check_axioms(build; expected=["SeedPrototype.h$(n)_ritz_witness"
        for (n, _) in WITNESS_SIZES])
    println("lake build (clean) passed in $(round(t; digits=1)) s; standard axioms only for ",
        join(names, ", "))

    full || return
    audit = lean_output(`$LAKE env lean AxiomAudit.lean`)
    names = check_axioms(audit; expected=["SeedPrototype.radixPack_extract",
        "SeedPrototype.weyl_counting_one", "SeedPrototype.h10_ritz_witness"])
    println("axiom audit: $(length(names)) declarations, standard axioms only")
    for file in sort(readdir(joinpath(SRC, "Bench")))
        t = @elapsed lean_output(`$LAKE env lean $(joinpath("Bench", file))`)
        println("Bench/$file checked in $(round(t; digits=1)) s")
    end
end

main("--full" in ARGS)
