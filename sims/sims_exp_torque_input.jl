using DifferentialEquations
using Plots
using Interpolations, Statistics
using MAT, YAML, JSON, HDF5

include("lib_types.jl")
include("lib_sim.jl")
include("lib_in_func.jl")
include("utils.jl")

Δt = 0.0001
default(lw=1.5)
# ----------------------------------------------------------------------
# Load experimenta points

fname = "exp_highest_cd61.h5"
fid = h5open("./sims/data-sim/" * fname, "r")
dst = "cd61_2_f"
dst = "cd61_3_f"
dst = "cd61_high_f"
dst = "cd61_high"

μ = mu_func(mean(fid[dst]["Watercut"][:]) / 100, mean(fid[dst]["T_3"][:]))
ρ = mean(fid[dst]["Density_oil"][:])
τ = fid[dst]["ESP_torque"][:]
x = fid[dst]["x"][:]
ω_ref = fid[dst]["ESP_rotation"][:]
P9_ref = fid[dst]["P_9"][:] * 100000 / ConvConst().H
P1_ref = fid[dst]["P_1"][:] * 100000 / ConvConst().H
DP_ref = P9_ref - P1_ref
q_ref = fid[dst]["Q_oil"][:]
a = fid[dst]["valve"][]
ω_t = fid[dst]["booster_oil"][]

# Generate a function for the torque
τf = linear_interpolation(x, τ)

# ----------------------------------------------------------------------
# Definitions

# Fluids
# Experimental mean of cd_i61 steady state files
oil = Fluid(1290, ρ, μ)

# Shaft
shaft = Shaft(0.2, 28.5, 1, 7850)
fit_mdl = ["e-2300-3300", "euler_fric"]
fit_mdl = ["e-2300-3300", "pinn_euler_fric"]
shaft = Shaft(28.5, 1, 7850, "./sims/models-exp-fit/torque_fit.yml", fit_mdl)

# Valve
fit_mdl = "cd61"
valve = Valve("./sims/models-exp-fit/valve_params_fit.yml", 1.0, fit_mdl)

# Twin-screw pump
fit_mdl = ["l-31.5_d-0.0762_eps-4.6e-05_h-0.5", "cd_i61"]
fit_mdl = ["l-31.5_d-0.0762_eps-4.6e-05_h-0.5", "pinn_cd_i61"]
fit_path = "./sims/models-exp-fit/twin_screw_fit.yml"
twin = Twin(ω_t, fit_path, fit_mdl)
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
fit_mdl = ["e-2300-3300", "pinn_euler_fric_loc"]
fit_path = "./sims/models-exp-fit/pump_fit.yml"
# p100_fit = Impeller(fit_path, fit_mdl, (108 + 60.5) / 2, (13 + 24) / 2, (108 - 60.5), oil, 77623.38803123798)
p100_fit = Impeller(fit_path, fit_mdl, (108 + 60.5) / 2, (13 + 24) / 2, (108 - 60.5), oil, 0.0)

# Systems
sys = System(oil, [pipe_up_twin, pipe_up, pipe_up, pipe_down, pipe_down_valve, pipe_down], shaft, p100_fit)
sys1 = System(oil, [pipe1_up, pipe1_down], shaft, p100_fit)

# Inputs
function steadystate(t, ω, A, A1, A2, T1, T2)
    if t <= 10
        return smooth_startup_10(t, τ[1])
    else
        return τf(t - 10)
    end
end

# inpt_w1 = Inputs([0.2, 30, 0.0, 40, 0, 10.0], in_sin3)
inpt_w1 = Inputs([0.2, 30, 0.0, 40, 0, 10.0], steadystate)

# ----------------------------------------------------------------------
# hqp1i1tcs_oil_p100_w1_turb
name = "hqvp3i1tcvfit"
tspan = (0.0, 10 + x[end])

u0 = zeros(14)

# With shock loss
sim = Sim(sys, u0, tspan, hqvp3i1tcvfit!, inpt_w1)
prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, AutoVern7(Rodas4()), reltol=1e-8, abstol=1e-8)

# Resample
t_0 = 0.0
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

plot(x .+ 10, ω_ref, label="Exp.")
plot!(t, ω * 60 / (2 * π) .- 60, label="Sim.")

plot(x .+ 10, P1_ref)
plot!(t, P[3, :] .+ 25)

plot(x .+ 10, P9_ref)
plot!(t, P[4, :] .+ 3.8)

plot(x .+ 10, q_ref / ρ)
plot!(t, Q[1, :])

# ----------------------------------------------------------------------
# hqp1i1tcs_oil_p100_w1_turb with real input
name = "hqvp1i1tcvfit"
tspan = (0.0, 10 + x[end])

u0 = zeros(6)

# With shock loss
sim1 = Sim(sys1, u0, tspan, hqvp1i1tcvfit!, inpt_w1)
prob = ODEProblem(sim1.model, sim1.u0, sim1.tspan, sim1)
# sol = solve(prob, AutoVern7(Rodas4()), reltol=1e-8, abstol=1e-8)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

# Resample
t_0 = 0.0
t_0 = 10.0
t = collect(t_0:Δt:sim1.tspan[2])
# t_0 = 8
# t = collect(t_0:Δt:10)
nsol = sol(t)

sim1.sys.pump.K

plot(x .+ 10, ω_ref, label="Exp.")
plot!(t, nsol[2, :] * 60 / (2 * π), label="Sim.")

plot(x .+ 10, P1_ref)
plot!(t, nsol[5, :], lw=2)

plot(x .+ 10, P9_ref)
plot!(t, nsol[6, :] .+ 3.9)

# plot(x .+ 10, q_ref / ρ)
plot(x .+ 10, q_ref)
plot!(t, nsol[1, :])

store_sim(name, sol, sim1, t, 10.0)

# ----------------------------------------------------------------------
# hqp1i1tcs_oil_p100_w1_turb with real input
name = "hqvp1i1lcvfit_p100_wocf"
tspan = (0.0, 10 + x[end])

u0 = zeros(6)

# With shock loss
sim1 = Sim(sys1, u0, tspan, hqvp1i1lcvfit!, inpt_w1)
prob = ODEProblem(sim1.model, sim1.u0, sim1.tspan, sim1)
# sol = solve(prob, AutoVern7(Rodas4()), reltol=1e-8, abstol=1e-8)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

# Resample
t_0 = 0.0
t_0 = 10.0
t = collect(t_0:Δt:sim1.tspan[2])
# t_0 = 8
# t = collect(t_0:Δt:10)
nsol = sol(t)

plot(x .+ 10, ω_ref, label="Exp.")
plot!(t, nsol[2, :] * 60 / (2 * π), label="Sim.")

plot(x .+ 10, P1_ref, label="Exp.")
plot!(t, nsol[5, :], lw=2, label="Sim.")

plot(x .+ 10, P9_ref, label="Exp.")
plot!(t, nsol[6, :] .+ 4.2, label="Sim.")

plot(x .+ 10, q_ref / ρ, label="Exp.")
# plot(x .+ 10, q_ref, label="Exp.")
plot!(t, nsol[1, :], label="Sim.")

store_sim(name, sol, sim1, t, 10.0)


# plot(x .+ 10, ω_ref, label="Exp.")
# plot!(t, nsol[2, :] * 60 / (2 * π) .- 60, label="Sim.")

# plot(x .+ 10, P1_ref)
# plot!(t, nsol[5, :] .+ 25, lw=2)

# plot(x .+ 10, P9_ref)
# plot!(t, nsol[6, :] .+ 3.8)

# # plot(x .+ 10, q_ref / ρ)
# plot(x .+ 10, q_ref)
# plot!(t, nsol[1, :] .- 1e-5)
