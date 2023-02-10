using DifferentialEquations
using Plots
using MAT, YAML, JSON, HDF5

include("lib_types.jl")
include("lib_sim.jl")
include("lib_in_func.jl")
include("utils.jl")

Δt = 0.0001

# ----------------------------------------------------------------------
# Load experimental points
exp_data = JSON.parsefile("./sims/exp-points/cd61.json")

exp_idx = "12"
μ = exp_data["mu"][exp_idx]
ρ = exp_data["Density_oil"][exp_idx]
τ = exp_data["ESP_torque"][exp_idx]
ω_ref = exp_data["ESP_rotation"][exp_idx]
P9_ref = exp_data["P_9"][exp_idx] * 100000 / ConvConst().H
P1_ref = exp_data["P_1"][exp_idx] * 100000 / ConvConst().H
DP_ref = P9_ref - P1_ref
q_ref = exp_data["Q_oil"][exp_idx]
a = exp_data["valve"][exp_idx]

# ----------------------------------------------------------------------
# Definitions

# Fluids
# Experimental mean of cd_i61 steady state files
oil = Fluid(1290, ρ, μ)

# Shaft
shaft = Shaft(0.2, 28.5, 1, 7850)
fit_mdl = ["e-2300-3300", "euler_fric"]
shaft = Shaft(28.5, 1, 7850, "./sims/models-exp-fit/torque_fit.yml", fit_mdl)

# Valve
fit_mdl = "cd61"
valve = Valve("./sims/models-exp-fit/valve_params_fit.yml", 1.0, fit_mdl)

# Twin-screw pump
fit_mdl = ["l-31.5_d-0.0762_eps-4.6e-05_h-0.5", "cd_i61"]
fit_path = "./sims/models-exp-fit/twin_screw_fit.yml"
twin = Twin(1280.0, fit_path, fit_mdl)
k_up_f = YAML.load_file(fit_path)[fit_mdl[1]][fit_mdl[2]]["klp"]

# Pipes
# Mauricio: 1.128+0.67+1.645+7.62 + 16+2.6
# Natan: Upstream: 31.5 downstream: 28
d_pipe = 3 * 25.4         # (mm) Pipe diameter
ϵ_pipe = 4.6e-05          # (m) Pipe rugosity

l_up, k_up = 31.5 / 3, k_up_f / 3
l_dw, k_dw = 28 / 3, 0.0 / 3

pipe_up = Pipe(d_pipe, l_up, ϵ_pipe, oil, k_up)
pipe_up_twin = Pipe(d_pipe, l_up, ϵ_pipe, oil, twin, k_up)
pipe_down_valve = Pipe(d_pipe, l_dw, ϵ_pipe, oil, valve)
pipe_down = Pipe(d_pipe, l_dw, ϵ_pipe, oil)

pipe1_up = Pipe(d_pipe, 31.5, ϵ_pipe, oil, twin, k_up_f)
pipe1_down = Pipe(d_pipe, 28.0, ϵ_pipe, oil, valve)

# Impellers
fit_mdl = ["e-2300-3300", "euler_fric_loc"]
fit_path = "./sims/models-exp-fit/pump_fit.yml"
p100_fit = Impeller(fit_path, fit_mdl, (108 + 60.5) / 2, (13 + 24) / 2, (108 - 60.5), oil, 77623.38803123798)

# Systems
sys = System(oil, [pipe_up_twin, pipe_up, pipe_up, pipe_down, pipe_down_valve, pipe_down], shaft, p100_fit)
sys1 = System(oil, [pipe1_up, pipe1_down], shaft, p100_fit)

# Inputs
function steadystate(t, ω, A, A1, A2, T1, T2)
    return smooth_startup_10(t, τ)
end

# inpt_w1 = Inputs([0.2, 30, 0.0, 40, 0, 10.0], in_sin3)
inpt_w1 = Inputs([0.2, 30, 0.0, 40, 0, 10.0], steadystate)

# ----------------------------------------------------------------------
# hqp1i1tcs_oil_p100_w1_turb
name = "hqvp3i1tcvfit"
tspan = (0.0, 15.0)

u0 = zeros(14)

# With shock loss
sim = Sim(sys, u0, tspan, hqvp3i1tcvfit!, inpt_w1)
prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, AutoVern7(Rodas4()), reltol=1e-8, abstol=1e-8)

# Resample
t_0 = 14.0
t_0 = 10.0
t = collect(t_0:Δt:sim.tspan[2])
# t_0 = 8
# t = collect(t_0:Δt:10)
nsol = sol(t)

Q_p, ω = view(nsol, 1, :), view(nsol, 2, :)
Q, P = view(nsol, 3:8, :), view(nsol, 9:14, :)

ω_ref
ω[end] * 60 / (2 * π)
q_ref
Q[end, end]
DP_ref
P[4, end] - P[3, end]
P1_ref
P[3, end]
P9_ref
P[4, end]


# ----------------------------------------------------------------------
# hqp1i1tcs_oil_p100_w1_turb
name = "hqvp1i1tcvfit"
tspan = (0.0, 15.0)

u0 = zeros(6)

# With shock loss
sim1 = Sim(sys1, u0, tspan, hqvp1i1tcvfit!, inpt_w1)
prob = ODEProblem(sim1.model, sim1.u0, sim1.tspan, sim1)
sol = solve(prob, AutoVern7(Rodas4()), reltol=1e-8, abstol=1e-8)

# Resample
t_0 = 14.0
t_0 = 10.0
t = collect(t_0:Δt:sim.tspan[2])
# t_0 = 8
# t = collect(t_0:Δt:10)
nsol = sol(t)

Q_p, ω = view(nsol, 1, :), view(nsol, 2, :)
Q, P = view(nsol, 3:4, :), view(nsol, 5:6, :)

ω_ref
ω[end] * 60 / (2 * π)
q_ref
Q[end, end]
DP_ref
P[2, end] - P[1, end]
P1_ref
P[1, end]
P9_ref
P[2, end]

# minimum([exp_data["booster_oil"][k] for k in keys(exp_data["booster_oil"])])

t = collect(10:Δt:sim.tspan[2])
store_sim(name, sol, sim, t, 10.0)

sim_dct = sim2dict(sim)

# Store info
f = open("tmp.yml", "w")
YAML.write(f, Dict(Symbol(name) => sim_dct))
close(f)
