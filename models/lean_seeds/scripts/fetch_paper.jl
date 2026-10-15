# Download the paper sources cited by the spin-1 studies (METHOD.md) into the gitignored
# upstream/, at the pinned commit, and check the hashes in results/spin1_kernel/
# paper-source-sha256.txt. Needed by scripts/seed_design/summarize.jl; network access required.
#
#     julia --project=. scripts/fetch_paper.jl
using Downloads, SHA
const MODEL = dirname(@__DIR__)
const UPSTREAM = joinpath(MODEL, "upstream")
const COMMIT="adc7f1241b42e322a6451854ab7e4b4c146bf78a"
const BASE="https://raw.githubusercontent.com/openai/math/$COMMIT/preprints/" *
    "The-periodic-spin-one-Haldane-gap-September-24-2026"
mkpath(UPSTREAM)
for name in ("initialization.tex","transfer.tex","thermal72.tex","trials.tex","bootstrap.tex")
    path=joinpath(UPSTREAM,name)
    Downloads.download("$BASE/build/$name",path)
    println(name," ",bytes2hex(sha256(read(path))))
end
for line in readlines(joinpath(MODEL, "results", "spin1_kernel", "paper-source-sha256.txt"))
    hash, path = split(line)
    bytes2hex(sha256(read(joinpath(MODEL, path)))) == hash || error("hash mismatch for $path")
end
println("paper sources match results/spin1_kernel/paper-source-sha256.txt")
