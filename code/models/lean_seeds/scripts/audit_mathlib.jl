# spin1_kernel: source-level Mathlib API audit (textual search with ripgrep, not a Lean proof of
# absence). Writes results/spin1_kernel/mathlib-audit.txt and mathlib-negative-search.txt.
# MATHLIB_DIR defaults to the Mathlib checkout of the repository's lean/ project, which pins the
# same toolchain (v4.35.0-rc3) and Mathlib commit c55e6e78. Needs `rg` and `git` on PATH.
#
#     julia scripts/audit_mathlib.jl          (from models/lean_seeds)
#
# This driver sits in scripts/, not scripts/spin1_kernel/: every .jl file under
# scripts/spin1_kernel/ is part of the source-hash set that check.jl compares with
# results/spin1_kernel/witnesses.json, so adding one there would invalidate that record.

const MODEL = dirname(@__DIR__)
const REPO = normpath(joinpath(MODEL, "..", ".."))
const PKG = get(ENV, "MATHLIB_DIR", joinpath(REPO, "lean", ".lake", "packages", "mathlib"))
const M = joinpath(PKG, "Mathlib")
const OUT = joinpath(MODEL, "results", "spin1_kernel")

isdir(M) || error("no Mathlib source at $M; set MATHLIB_DIR")

"""
    rg(args; allow_empty = false) -> Vector{String}

Lines of `rg -n --no-heading args` with the Mathlib prefix removed. A search without hits
(exit status 1) is an error unless `allow_empty`.
"""
function rg(args; allow_empty = false)
    io = IOBuffer()
    process = run(pipeline(ignorestatus(`rg -n --no-heading $args`); stdout = io))
    status = process.exitcode
    status == 0 || (status == 1 && allow_empty) || error("rg $args failed with status $status")
    lines = split(String(take!(io)), '\n'; keepempty = false)
    return [replace(line, M * "/" => "") for line in lines]
end

file(path) = joinpath(M, path)

audit = String[strip(read(`git -C $PKG rev-parse HEAD`, String))]
append!(audit, rg(["-e", "def eigenvalues|eigenvalues₀_antitone|spectral_theorem",
    "-e", "mulVec_eigenvectorBasis|trace_eq_sum_eigenvalues|eigenvalues_eq_eigenvalues_iff",
    file("Analysis/Matrix/Spectrum.lean")]))
append!(audit, rg(["-e", "hasEigenvalue_iInf|hasEigenvalue_iSup|rayleighQuotient_le_norm",
    "-e", "norm_eq_iSup_rayleighQuotient|hasEigenvector_of_isMinOn",
    file("Analysis/InnerProductSpace/Rayleigh.lean")]))
append!(audit, rg(["posSemidef_iff_eigenvalues_nonneg|eigenvalues_nonneg",
    file("Analysis/Matrix/PosDef.lean")]))
append!(audit, rg(["posSemidef_iff_dotProduct_mulVec|dotProduct_mulVec_nonneg|trace_nonneg",
    file("LinearAlgebra/Matrix/PosDef.lean")]))
append!(audit, first(rg(["map_sum_le", file("Analysis/Convex/Jensen.lean")]), 8))
append!(audit, rg(["convexOn_exp", file("Analysis/Convex/SpecificFunctions/Basic.lean")]))
append!(audit, rg(["tracePositiveLinearMap|lemma le_iff", file("Analysis/Matrix/Order.lean")]))
write(joinpath(OUT, "mathlib-audit.txt"), join(audit, '\n') * '\n')

# Negative search for min-max / Weyl / interlacing; hits are listed for manual triage.
negative = rg(["-i", "courant|fischer|interlac|peierls|weyl.{0,50}inequal|min.max.{0,50}eigen",
    M, "--glob", "*.lean", "--glob", "!**/RootSystem/**", "--glob", "!**/Lie/**",
    "--glob", "!**/Cartan/**"]; allow_empty = true)
write(joinpath(OUT, "mathlib-negative-search.txt"), isempty(negative) ? "" :
    join(negative, '\n') * '\n')
println("wrote mathlib-audit.txt ($(length(audit)) lines) and mathlib-negative-search.txt ",
    "($(length(negative)) hits)")
