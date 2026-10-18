using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Blocks as Blocks
using OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using Test

# Helpers shared by more than one test file.
simulate(model, T = 1.0) =
    solve(ODEProblem(mtkcompile(model), [], (0.0, T)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
const T1 = range(0.0, 1.0, length = 11)

solve_tp(model, T = 1.0) =
    solve(ODEProblem(mtkcompile(model), [], (0.0, T)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
const TT = range(0.0, 1.0, length = 11)

# A modulus read from a signal is an unknown, which MTK won't divide by, so the model
# keeps an algebraic equation: a DAE that needs a mass-matrix solver, and a nonlinear
# solve for its initial values when the loop is nonlinear.
solve_dae(model) =
    solve(ODEProblem(mtkcompile(model), [], (0.0, 1.0)), Rodas5P(); abstol = 1e-10, reltol = 1e-10)

# Se(E) ⇀ 1 ⇀ two-port ⇀ load, as one model.
function driven(twoport, load; E)
    @named se = Se(E)
    @named src = Junction1(se => :pout, :out)
    System([connect(src.out, twoport.p1), connect(twoport.p2, load.power)], t;
           name = :m, systems = [src, twoport, load])
end

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
    @testset "Graph" begin
        include("test_graph.jl")
    end
    @testset "Regression" begin
        include("regression/helpers.jl")
        for f in filter(startswith("test_"), readdir(joinpath(@__DIR__, "regression")))
            include(joinpath("regression", f))
        end
    end
end
