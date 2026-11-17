# Coupled-seed points certified by certify.jl and checked by check.jl (included by
# src/RitzSeeds.jl): transverse field g and seed (n₀ = RITZ_N0 rungs, β₀). β₀ is the best point
# of the floating-point scan of scripts/ritz_scan.jl on a 1/32 grid (see METHOD.md).
const RITZ_POINTS = [(; g = 3 // 1, n0 = RITZ_N0, β0 = 29 // 16),
    (; g = 4 // 1, n0 = RITZ_N0, β0 = 11 // 8), (; g = 11 // 4, n0 = RITZ_N0, β0 = 69 // 32)]
