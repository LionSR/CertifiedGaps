# Download the paper sources cited by the spin-1 studies (METHOD.md) into the gitignored
# upstream/ and check the hashes in results/spin1_kernel/paper-source-sha256.txt. Needed by
# scripts/seed_design/summarize.jl; network access required. The fetcher itself is the shared
# `tools/fetch_upstream.jl --paper` (one place for the pinned commit); this name is kept because
# the hash-pinned `src/Common.jl` documents `upstream/` by it.
#
#     julia scripts/fetch_paper.jl        (equivalently, from the repository root:
#                                          julia tools/fetch_upstream.jl --paper)
run(`$(Base.julia_cmd()) --startup-file=no $(joinpath(@__DIR__, "..", "..", "..", "tools",
    "fetch_upstream.jl")) --paper $ARGS`)
