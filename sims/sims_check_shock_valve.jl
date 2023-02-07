using DifferentialEquations
using Plots
using MAT, YAML

include("lib_types.jl")
include("lib_sim.jl")
include("lib_in_func.jl")
include("utils.jl")


Δt = 0.0001

# ----------------------------------------------------------------------
# Definitions

# Fluids
oil = Fluid(1290, 805, 50e-3)
oil = Fluid(1290, 805, 150e-3)

# Shaft
shaft = Shaft(0.1, 28.5, 1, 7850)

# Valve
valve = Valve("./sims/cfg/valve_params.yml", 1.0)

# Pipes
pipe_oil = Pipe(3 * 25.4, 7.5, 0.1 * 1e-3, oil)
pipe_oil_valve = Pipe(3 * 25.4, 7.5, 0.1 * 1e-3, oil, valve)

# Impellers
p100_oil = Impeller(0.002, -1.7317, 108, 0.2, (108 + 60.5) / 2, (13 + 24) / 2, (108 - 60.5), 100, 0.032, oil)
# Systems
sys = System(oil, pipe_oil, shaft, p100_oil)
sys_valve = System(oil, pipe_oil_valve, shaft, p100_oil)

# Inputs
inpt_w3 = Inputs([3.0, 30, 1e-3, 40, 0.25, 2.0], in_sin2)
inpt_w1 = Inputs([1.0, 30, 0.0, 40, 0, 10.0], in_sin2)

# ----------------------------------------------------------------------
# hqp1i1tcs_oil_p100_w1_turb
name = "hqp1i1tcs_oil_p100_w1_turb"
tspan = (0.0, 15.0)

#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# With shock loss
sim = Sim(sys, u0, tspan, hqp1i1tcs!, inpt_w1)
prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sols = solve(prob, reltol=1e-9, abstol=1e-9)

# Without shock loss
sim = Sim(sys, u0, tspan, hqp1i1tc!, inpt_w1)
prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-9, abstol=1e-9)

# Resample
t_0 = 10.0
t = collect(t_0:Δt:sim.tspan[2])
nsols = sols(t)
nsol = sol(t)

# Flowrate
plot(t, nsols[3, :], label="W shock")
plot!(t, nsol[3, :], label="Wo shock")
# Pressure
plot(t, nsols[5, :], label="W shock")
plot!(t, nsol[5, :], label="Wo shock")
# Reynolds
plot(sols.t, ReQ.(sim.sys.fluid.ρ, sols[3, :], sim.sys.pipe.d, sim.sys.fluid.μ))

# ----------------------------------------------------------------------
# hqvp1i1tcs_oil_p100_w1_turb_valve
name = "hqvp1i1tcs_oil_p100_w1_turb_valve"
tspan = (0.0, 15.0)

#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# With shock loss
sim = Sim(sys_valve, u0, tspan, hqvp1i1tcs!, inpt_w1)
prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
solv = solve(prob, reltol=1e-9, abstol=1e-9)
# For the valve pressure, I also checked with the tool below and I got the same 
# result. https://www.swagelok.com/en/toolbox/cv-calculator

plot(solv)
nsolv = solv(t)

# Flowrate
plot(t, nsolv[3, :], label="W valve")
plot!(t, nsols[3, :], label="W shock")
plot!(t, nsol[3, :], label="Wo shock")
# Pressure
plot(t, nsolv[5, :], label="W valve")
plot!(t, nsols[5, :], label="W shock")
plot!(t, nsol[5, :], label="Wo shock")

plot(t, nsolv[6, :], label="W valve")
plot!(t, nsols[6, :], label="W shock")
plot!(t, nsol[6, :], label="Wo shock")
# Reynolds
plot(sols.t, ReQ.(sim.sys.fluid.ρ, sols[3, :], sim.sys.pipe.d, sim.sys.fluid.μ))

