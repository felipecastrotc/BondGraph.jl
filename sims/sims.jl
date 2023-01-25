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
oil = Fluid(1290, 850, 800e-3)

# Shaft
shaft = Shaft(0.1, 28.5, 1, 7850)
shaft_b = Shaft(0.1, 100, 2, 7850)

# Pipes
pipe_oil = Pipe(3 * 25.4, 7.5, 0.1 * 1e-3, oil)

# Impellers
p100_oil = Impeller(0.002, -1.7317, 108, 0.2, (108 + 60.5) / 2, (13 + 24) / 2, (108 - 60.5), oil)
# Systems
sys = System(oil, pipe_oil, shaft, p100_oil)
sys_bshat = System(oil, pipe_oil, shaft_b, p100_oil)

# Inputs
inpt_w1 = Inputs([1.0, 30, 1e-3, 40, 0.25, 1.0], in_sin)
inpt_w3 = Inputs([3.0, 30, 1e-3, 40, 0.25, 2.0], in_sin2)
inpt_w1 = Inputs([1.0, 30, 1e-4, 40, 0.25, 2.0], in_sin2)

# ----------------------------------------------------------------------
# hqp1i1_oil_p100_bshaft_w0
name = "hqp1i1_oil_p100_bshaft_w0"
# the bshaft means it uses a larger shaft for better results on simulation
tspan = (0.0, 6.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys_bshat, u0, tspan, hqp1i1!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

t_0 = 2.0
t = collect(t_0:Δt:sim.tspan[2])
store_sim(name, sol, sim, t, t_0)

plot(sol)
# ----------------------------------------------------------------------
# hqp1i1_oil_p100_w0
name = "hqp1i1_oil_p100_w0"
tspan = (0.0, 6.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys, u0, tspan, hqp1i1!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-9, abstol=1e-9)

t_0 = 2.0
t = collect(t_0:Δt:sim.tspan[2])
store_sim(name, sol, sim, t, t_0)

# ----------------------------------------------------------------------
# hqp1i1_oil_p100_w1
name = "hqp1i1_oil_p100_w1"
tspan = (0.0, 6.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys, u0, tspan, hqp1i1!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

t_0 = 2.0
t = collect(t_0:Δt:sim.tspan[2])
store_sim(name, sol, sim, t, t_0)

# ----------------------------------------------------------------------
# hqp1i1_oil_p100_w3
name = "hqp1i1_oil_p100_w3"
tspan = (0.0, 5.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys, u0, tspan, hqp1i1!, inpt_w3)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

t_0 = 2.0
t = collect(t_0:Δt:sim.tspan[2])
store_sim(name, sol, sim, t, t_0)

plot(sol)