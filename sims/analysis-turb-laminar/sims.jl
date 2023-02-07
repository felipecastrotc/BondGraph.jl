using DifferentialEquations
using Plots
using MAT, YAML

include("../lib_types.jl")
include("../lib_sim.jl")
include("../lib_in_func.jl")
include("../utils.jl")

function gen_plots(sim, sol, idx_start=1000, visc="", save=false)
    path = "./sims/analysis-turb-laminar/figs/"

    # Get sys info
    s = sim.sys
    μ, ρ = s.fluid.μ, s.fluid.ρ
    l, d, ϵ, I = s.pipe.l, s.pipe.d, s.pipe.ϵ, s.pipe.I

    # Resample system
    ts = 0:0.001:sim.tspan[2]
    nsol = sol(ts)
    Q = nsol[1, :]
    # Calculate reynolds over time
    Re = ReQ.(ρ, Q, d, μ)

    # Plot input
    display(plot(ts, inpt_w1.func.(ts, inpt_w1.p...), title=visc, ylabel="Torque (N*m)", xlabel="Time (s)"))
    if save
        savefig(path * "input-func_" * visc * ".svg")
    end
    # Plot Reynolds
    display(plot(ts, Re, ylabel="Reynolds (-)", xlabel="Time (s)", title=visc))
    if save
        savefig(path * "reynolds_" * visc * ".svg")
    end

    # Friction factor along the time
    slc = idx_start:length(Re)
    plot(ts[slc], 64.0 ./ Re[slc], label="laminar")
    plot!(ts[slc], cheng.(abs.(Re), ϵ, d)[slc], label="cheng")
    plot!(ts[slc], haaland.(Q, d, ϵ, μ, ρ)[slc], label="haaland")
    xlabel!("Time (s)")
    title!(visc)
    display(ylabel!("Friction factor (-)"))
    if save
        savefig(path * "friction_" * visc * ".svg")
    end

    plot(ts[slc], cheng.(abs.(Re), ϵ, d)[slc], label="cheng")
    plot!(ts[slc], haaland.(Q, d, ϵ, μ, ρ)[slc], label="haaland")
    xlabel!("Time (s)")
    title!(visc)
    display(ylabel!("Friction factor (-)"))
    if save
        savefig(path * "friction_turb_" * visc * ".svg")
    end

    plot(ts[slc], Q[slc], label="Q", title="Flow-rate along time")
    ylabel!("Q (m^3/s)")
    title!(visc)
    display(xlabel!("Time (s)"))
    if save
        savefig(path * "flow-rate_" * visc * ".svg")
    end

end

Δt = 0.0001

# ----------------------------------------------------------------------
# Definitions

# Fluids
hvisc = Fluid(1290, 850, 800e-3)
mvisc = Fluid(1290, 850, 60e-3)
lvisc = Fluid(1290, 850, 5e-3)
wvisc = Fluid(1290, 998, 1e-3)

# Shaft
shaft = Shaft(0.1, 28.5, 1, 7850)

# Pipes
pipe_hvisc = Pipe(3 * 25.4, 7.5, 0.1 * 1e-3, hvisc)
pipe_mvisc = Pipe(3 * 25.4, 7.5, 0.1 * 1e-3, mvisc)
pipe_lvisc = Pipe(3 * 25.4, 7.5, 0.1 * 1e-3, lvisc)
pipe_wvisc = Pipe(3 * 25.4, 7.5, 0.1 * 1e-3, wvisc)

# Impellers
p100_hvisc = Impeller(0.002, -1.7317, 108, 0.2, (108 + 60.5) / 2, (13 + 24) / 2, (108 - 60.5), hvisc)
p100_mvisc = Impeller(p100_hvisc, hvisc)
p100_lvisc = Impeller(p100_hvisc, lvisc)
p100_wvisc = Impeller(p100_hvisc, wvisc)

# Systems
sys_hvisc = System(hvisc, pipe_hvisc, shaft, p100_hvisc)
sys_mvisc = System(mvisc, pipe_mvisc, shaft, p100_mvisc)
sys_lvisc = System(lvisc, pipe_lvisc, shaft, p100_lvisc)
sys_wvisc = System(wvisc, pipe_lvisc, shaft, p100_wvisc)

# Inputs
inpt_w3 = Inputs([3.0, 30, 1e-3, 40, 0.25, 2.0], in_sin2)
inpt_w1 = Inputs([1.0, 30, 0.0, 40, 0, 10.0], in_sin2)

# ----------------------------------------------------------------------
# High viscosity
name = "high_viscosity"
tspan = (0.0, 20.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys_hvisc, u0, tspan, hqp1i1th!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)
gen_plots(sim, sol, 1000, name, true)

# ----------------------------------------------------------------------
# Medium viscosity

name = "medium_viscosity"
tspan = (0.0, 20.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys_mvisc, u0, tspan, hqp1i1th!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)
gen_plots(sim, sol, 1000, name, true)

# ----------------------------------------------------------------------
# Low viscosity
name = "low_viscosity"
tspan = (0.0, 20.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys_lvisc, u0, tspan, hqp1i1th!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)
gen_plots(sim, sol, 1000, name, false)
# gen_plots(sim, sol, 1000, name, true)

# ----------------------------------------------------------------------
# Low viscosity
name = "water"
tspan = (0.0, 20.0)
#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
sim = Sim(sys_wvisc, u0, tspan, hqp1i1th!, inpt_w1)

prob = ODEProblem(sim.model, sim.u0, sim.tspan, sim)
sol = solve(prob, reltol=1e-8, abstol=1e-8)
gen_plots(sim, sol, 1000, name, true)
