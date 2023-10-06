using BondGraph
using Test

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

@info "All tests completed!"

# TODO: higher level tests like in the example