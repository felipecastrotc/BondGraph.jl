using BondGraphToolkit
using OrdinaryDiffEqTsit5
using Test

@testset "BondGraphToolkit" begin
    @testset "Core" begin
        include("test_core.jl")
    end
    @testset "Mechanical" begin
        include("test_mechanical.jl")
    end
    @testset "Multi-junction" begin
        include("test_multi_junction.jl")
    end
    @testset "Two-ports" begin
        include("test_two_ports.jl")
    end
    @testset "Signals" begin
        include("test_signals.jl")
    end
    @testset "Adapters" begin
        include("test_adapters.jl")
    end
    @testset "Equations" begin
        include("test_equations.jl")
    end
    @testset "Regression" begin
        include("regression/helpers.jl")
        for f in filter(startswith("test_"), readdir(joinpath(@__DIR__, "regression")))
            include(joinpath("regression", f))
        end
    end
end
