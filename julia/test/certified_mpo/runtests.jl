# CertifiedMPO tests, ported from spin2/test. Each file runs in its own module, as it would
# standalone, so that test-local helpers of the same name do not collide.
for name in (
    "mpo_norms",
    "mpo_squaring",
    "mpo_products",
    "charge_blocks",
    "mpo_operator_error",
    "mpo_gauge_search",
    "adjoint_trace",
    "schatten_bridge",
    "positive_bounds",
    "polynomial_bounds",
    "polynomial_mpo",
)
    path = joinpath(@__DIR__, name * ".jl")
    Core.eval(Main, :(module $(Symbol("CertifiedMPOTest_", name))
    include($path)
    end))
end
