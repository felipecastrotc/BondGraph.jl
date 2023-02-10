# =====================================================================
# This file was used to generate the intial results on the PINN dataset
# It is being kept for reproducibility
# =====================================================================

using DifferentialEquations
using Plots
using Statistics, Interpolations
using MAT, YAML, HDF5

include("lib_types.jl")
include("lib_sim.jl")
include("lib_in_func.jl")
include("utils.jl")

Δt = 0.0001

# ----------------------------------------------------------------------
# Load experimental data

fname = "exp_highest_cd61.h5"
fid = h5open("./sims/data-sim/" * fname, "r")
dst = "cd61_high_f"

# Load properties
μ = mu_func(mean(fid[dst]["Watercut"][:]) / 100, mean(fid[dst]["T_3"][:]))
ρ = mean(fid[dst]["Density_oil"][:])
ω_t = fid[dst]["booster_oil"][]
τ = fid[dst]["ESP_torque"][:]
x = fid[dst]["x"][:]
fid.close()

# Store the input info for reproducibility
fname = "reproducibility_torque_sig_highest_cd61.h5"
fid = h5open("./sims/data-sim/"*fname, "w")
fid["ESP_torque"] = τ

# μ = 0.14293630582220418      # Kept for reproducibility
# ρ = 882.1995559385915        # Kept for reproducibility
# ω_t = 1280.0                 # Kept for reproducibility

# Generate a function for the torque
τf = linear_interpolation(x, τ)

# ----------------------------------------------------------------------
# Definitions

# Fluids
# Experimental mean of cd_i61 steady state files
oil = Fluid(1290, ρ, μ)

# Shaft
fit_mdl = ["e-2300-3300", "euler_fric"]
shaft = Shaft(28.5, 1, 7850, "./sims/models-exp-fit/torque_fit.yml", fit_mdl)

# Valve
fit_mdl = "cd61"
valve = Valve("./sims/models-exp-fit/valve_params_fit.yml", 1.0, fit_mdl)

# Twin-screw pump
fit_mdl = ["l-31.5_d-0.0762_eps-4.6e-05_h-0.5", "cd_i61"]
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
pipe_up = Pipe(d_pipe, 31.5, ϵ_pipe, oil, twin, k_up_f)
pipe_down = Pipe(d_pipe, 28.0, ϵ_pipe, oil, valve)

# Impellers
fit_mdl = ["e-2300-3300", "euler_fric_loc"]
fit_path = "./sims/models-exp-fit/pump_fit.yml"
p100_fit = Impeller(fit_path, fit_mdl, (108 + 60.5) / 2, (13 + 24) / 2, (108 - 60.5), oil, 0.0)

# Systems
sys1 = System(oil, [pipe_up, pipe_down], shaft, p100_fit)

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
# hqp1i1_oil_p100_bshaft_w0
name = "hqp1i1_oil_p100_bshaft_w0"
# the bshaft means it uses a larger shaft for better results on simulation
tspan = (0.0, 5.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys_bshat, u0, tspan, hqp1i1!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

t = collect(1:Δt:2)
store_sim(name, sol, sim, t, 1.0)

# ----------------------------------------------------------------------
# hqp1i1_oil_p100_w0
name = "hqp1i1_oil_p100_w0"
tspan = (0.0, 5.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys, u0, tspan, hqp1i1!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

t = collect(1:Δt:2)
store_sim(name, sol, sim, t, 1.0)

# ----------------------------------------------------------------------
# hqp1i1_oil_p100_w1
name = "hqp1i1_oil_p100_w1"
tspan = (0.0, 5.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys, u0, tspan, hqp1i1!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

t = collect(1:Δt:sim.tspan[2])
store_sim(name, sol, sim, t, 1.0)

# ----------------------------------------------------------------------
# hqp1i1_oil_p100_w3
name = "hqp1i1_oil_p100_w3"
tspan = (0.0, 5.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys, u0, tspan, hqp1i1!, inpt_w3)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

t = collect(1:Δt:sim.tspan[2])
store_sim(name, sol, sim, t, 1.0)

plot(sol)