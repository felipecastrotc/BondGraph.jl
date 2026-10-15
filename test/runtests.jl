using BondGraph
using Test

# `Pkg.test(test_args = ["regenerate"])` rewrites the golden data and stops.
if "regenerate" in ARGS
    include("regression/generate_references.jl")
    exit()
end

# Test types and basic utilities
@testset "Test base" begin
    include("unit-tests/test_types.jl")
    include("unit-tests/test_graph.jl")
end

# Utilities for bond graphs and connections
@testset "Test utilities" begin
    include("unit-tests/test_bg_utils.jl")
    include("unit-tests/test_connection_utils.jl")
end

# System and port-related tests
@testset "Test ports and system " begin
    include("unit-tests/test_ports.jl")
    include("unit-tests/test_system_utils.jl")
end

# Macro-related tests
@testset "Test macros " begin
    include("unit-tests/test_macros.jl")
end

# Golden-data regression of every example (see test/regression/helpers.jl)
@testset "Regression" begin
    include("regression/helpers.jl")
    for f in filter(startswith("test_"), readdir(joinpath(@__DIR__, "regression")))
        include(joinpath("regression", f))
    end
end

@info "All tests completed!"
