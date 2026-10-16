using BondGraphToolkit
using OrdinaryDiffEqTsit5
using Test

# Examples not yet ported to v0.3, by golden file → ticket in .scratch/v0.3-impl/issues.
const PENDING = Dict(
    "test_mechanical_msd_1dof.jl" => "02",
    "test_mechanical_duffing.jl" => "02",
    "test_electrical_dc_motor.jl" => "04",
)

@testset "BondGraphToolkit" begin
    @testset "Core" begin
        include("test_core.jl")
    end
    @testset "Multi-junction" begin
        include("test_multi_junction.jl")
    end
    @testset "Regression" begin
        include("regression/helpers.jl")
        for f in filter(startswith("test_"), readdir(joinpath(@__DIR__, "regression")))
            haskey(PENDING, f) ? (@info "Skipping $f until ticket $(PENDING[f])") : include(joinpath("regression", f))
        end
    end
end
