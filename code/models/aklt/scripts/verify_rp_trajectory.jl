#!/usr/bin/env julia

# Standalone exact arithmetic checks. No model environment or external package is required.
# These verify an obstruction and conditional estimates, not an AKLT spectral-gap certificate.
include(joinpath(@__DIR__, "..", "test", "rigorous_checks.jl"))
