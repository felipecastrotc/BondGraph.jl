using DifferentialEquations
using Plots
using MAT, YAML

include("lib_types.jl")
include("lib_sim.jl")
include("lib_in_func.jl")
include("lib_flow.jl")
include("utils.jl")


1.128 + 0.67 + 1.645 + 7.62 + 16 + 2.6

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
# Progressing cavity pump
cavity = Twin(1.0, 2e-6, 0.05, 0.0)

1300*d = 0.05
0.05/1300
pcavity(0.04, 3e6, 0.05, 1.0, 1)/ConvConst().H

# Pipes
# Mauricio: 1.128+0.67+1.645+7.62 + 16+2.6
# Natan: Upstream: 31.5 downstream: 28
d_pipe = 3 * 25.4         # (mm) Pipe diameter
ϵ_pipe = 0.1 * 1e-3       # (m) Pipe rugosity
pipe_up = Pipe(d_pipe, 31.5, ϵ_pipe, oil)
pipe_up_cavity = Pipe(d_pipe, 31.5, ϵ_pipe, oil, cavity)
pipe_down = Pipe(d_pipe, 11.043, ϵ_pipe, oil)
pipe_down_valve = Pipe(d_pipe, 16.957, ϵ_pipe, oil, valve)

# Impellers
p100_oil = Impeller(0.002, -1.7317, 108, 0.2, (108 + 60.5) / 2, (13 + 24) / 2, (108 - 60.5), 100, 0.032, oil)
# Systems
sys = System(oil, [pipe_up, pipe_down, pipe_down_valve], shaft, p100_oil)
sys_cav = System(oil, [pipe_up_cavity, pipe_down, pipe_down_valve], shaft, p100_oil)

# Inputs
inpt_w3 = Inputs([3.0, 30, 1e-3, 40, 0.25, 2.0], in_sin2)
inpt_w1 = Inputs([1.0, 10, 0.0, 40, 0, 10.0], in_sin2)

# ----------------------------------------------------------------------
# hqp1i1tcs_oil_p100_w1_turb
name = "hqvp2i1tcs!_oil_p100_w1_turb"
tspan = (0.0, 15.0)

#      Qi,   ω,  Q1,  Q3,  Q5,  P2,  P4, P6
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# With shock loss
sim = Sim(sys, u0, tspan, hqvp2i1tcs!, inpt_w1)
prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

# Resample
t_0 = 10.0
t = collect(t_0:Δt:sim.tspan[2])
nsol = sol(t)

plot(t, nsol[1, :], label="W shock")
# Flowrate
plot(t, nsol[3, :], label="W shock")
plot!(t, nsol[4, :], label="W shock")
plot!(t, nsol[5, :], label="W shock")
# Pressure
plot(t, nsol[6, :], label="W shock")
plot!(t, nsol[7, :], label="Wo shock")
plot!(t, nsol[8, :], label="Wo shock")
plot(t, nsol[2, :], label="Wo shock")
# Reynolds
plot(sols.t, ReQ.(sim.sys.fluid.ρ, sol[3, :], sim.sys.pipe[1].d, sim.sys.fluid.μ))

# ----------------------------------------------------------------------
# hqvp1i1tcs_oil_p100_w1_turb_valve
name = "hqvp2i1tcsc_oil_p100_w1_turb_valve"
tspan = (0.0, 15.0)

#      Qi,   ω,  Q1,  Q3,  Q5,  P2,  P4, P6
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# With shock loss
sim = Sim(sys_cav, u0, tspan, hqvp2i1tcsc!, inpt_w1)
prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
solc = solve(prob, reltol=1e-9, abstol=1e-9)
# For the valve pressure, I also checked with the tool below and I got the same 
# result. https://www.swagelok.com/en/toolbox/cv-calculator

plot(solc)
t = collect(0:Δt:sim.tspan[2])
nsolc = solc(t)

plot(t, nsolc[2, :], label="W valve")
# Flowrate
plot(t, nsolc[3, :], label="W valve")
plot!(t, nsolc[4, :], label="W valve")
plot!(t, nsolc[5, :], label="W valve")
# Pressure
plot(t, nsolc[6, :], label="W valve")
plot!(t, nsolc[7, :], label="W shock")
plot!(t, nsolc[8, :], label="Wo shock")

P_0 = pcavity.(nsolc[3, :], oil.μ, cavity.ω, cavity.k_1, cavity.k_2, 0.0)
plot(P_0)


# Reynolds
plot(sols.t, ReQ.(sim.sys.fluid.ρ, sols[3, :], sim.sys.pipe.d, sim.sys.fluid.μ))
