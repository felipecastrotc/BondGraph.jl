# Golden-data regression helpers. The only place that knows the CSV layout and
# how states are named; adapt `state_names` if MTK renames states after an upgrade.
using BondGraph: t, D  # @oneport and the examples expect these in scope
using DelimitedFiles
using DifferentialEquations
using Test

const REF_DIR = joinpath(@__DIR__, "reference")
const TS = range(0.0, 10.0, length = 21)  # sample times shared by every case
const REGENERATE = Ref(false)             # set by generate_references.jl

state_names(prob) = string.(states(prob.f.sys))

function simulate(prob)
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    sol.retcode == ReturnCode.Success || error("solve failed: $(sol.retcode)")
    return state_names(prob), reduce(hcat, sol(TS).u)'  # rows = TS, cols = states
end

# `build()` returns (sys, ["case" => ODEProblem, ...]); one CSV per case.
function check_example(example, build)
    sys, cases = build()
    if REGENERATE[]
        write(joinpath(REF_DIR, "$(example)_equations.txt"), join(string.(equations(sys)), "\n") * "\n")
        for (case, prob) in cases
            names, data = simulate(prob)
            writedlm(joinpath(REF_DIR, "$case.csv"), [permutedims(["t"; names]); [TS data]], ',')
        end
        return
    end
    @testset "$case" for (case, prob) in cases
        names, data = simulate(prob)
        ref, header = readdlm(joinpath(REF_DIR, "$case.csv"), ',', Float64; header = true)
        header = vec(header)
        @test ref[:, 1] ≈ TS
        for (k, name) in enumerate(header)
            k == 1 && continue  # `t` column
            j = findfirst(==(name), names)
            if j === nothing
                @error "$case: reference state `$name` missing from simulated system (has $names)"
                @test false
                continue
            end
            ok = all(isapprox.(data[:, j], ref[:, k]; rtol = 1e-6, atol = 1e-9))
            ok || @error "$case: state `$name` deviates, max abs diff $(maximum(abs, data[:, j] - ref[:, k]))"
            @test ok
        end
        extra = setdiff(names, header)
        isempty(extra) || @warn "$case: states not in reference: $extra"
    end
end
