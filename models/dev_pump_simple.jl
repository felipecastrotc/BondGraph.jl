using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations
using ModelingToolkit

# ==========================================================
# Used alias
@variables ω(t), Q(t)

# ==========================================================
# New test with multiple outputs for pump leakage

# ------------------------------------------------------------------------------
# Elements
@parameters p, po, τ, Qₚ, ρ

@named m = Mass()       # General inertial element
@named d = Damper()     # General damping element
@named dv = Damper()     # General damping element
@named s = Spring()     # General compliance element
@named P = Se(p)        # Inlet pressure
@named Qp = Sf(Qₚ)        # Positive displacement pump
@named Po = Se(po)        # Inlet pressure
@named T = Se(τ)        # Shaft input torque

# Pipeline elements
@named plm = Mass(m = rpo["Iₚᵢ"])
@named pld = Damper(c = rpo["fₚᵢ"])     # General damping element
@named pls = Spring(k = -1 / rpo["kᵤ"])     # General compliance element

# ------------------------------------------------------------------------------
# Impeller Junctions

@named lek = Junction1([-1, d])     # Leakage setting

@named jus = Junction0()     # Coupling for leakage

@named imp = Junction1([-1, m], [-1, d])     # Impeller

@named olt = Junction1([-1, d])     # Shock losses
@named jds = Junction0()            # Coupling for leakage

# ------------------------------------------------------------------------------
# Axis Junctions
@named ax = Junction1(T, [-1, m], [-1, d])      # Axis parameters

# ------------------------------------------------------------------------------
# Pipe Junctions
@named us11 = Junction1([-1, plm], [-1, pld], P)     # Inlet pipe
@named us10 = Junction0([-1, pls])
# @named us10 = Junction0()
@named us21 = Junction1([-1, plm], [-1, pld])     # Inlet pipe
@named us20 = Junction0([-1, pls])
# @named us20 = Junction0()
@named us31 = Junction1([-1, plm], [-1, pld])     # Inlet pipe
@named us30 = Junction0([-1, pls])
@named us41 = Junction1([-1, plm], [-1, pld])     # Inlet pipe
@named us40 = Junction0([-1, pls])
# @named us30 = Junction0()
plist = [us11, us10, us21, us20, us31, us30, us41, us40];

@named ds10 = Junction0([-1, pls])
# @named ds10 = Junction0()
@named ds11 = Junction1([-1, plm], [-1, pld])     # Outlet pipe
@named ds20 = Junction0([-1, pls])
# @named ds20 = Junction0()
@named ds21 = Junction1([-1, plm], [-1, pld])     # Outlet pipe
@named ds30 = Junction0([-1, pls])
# @named ds30 = Junction0()
@named ds31 = Junction1([-1, plm], [-1, pld])     # Outlet pipe
@named ds40 = Junction0([-1, pls])
# @named ds30 = Junction0()
@named ds41 = Junction1([-1, plm], [-1, pld], Po)     # Outlet pipe
push!(plist, ds11, ds10, ds21, ds20, ds31, ds30, ds41, ds40);

# Tank junction
@named tnkin = Junction0(P)
@named tnkout = Junction0(Po)

# ------------------------------------------------------------------------------
# MGY
@parameters γa, γb

gyp = ρ * (γa * ω - γb * Q)
@named agy = mGY(g = gyp)

# ------------------------------------------------------------------------------
# Build connections

# Pipeline assembly
cons = [
    connect(us11.power, us10.power),
    connect(us10.power, us21.power),
    connect(us21.power, us20.power),
    connect(us20.power, us31.power),
    connect(us31.power, us30.power),
    connect(us30.power, us41.power),
    connect(us41.power, us40.power),
]
push!(
    cons,
    connect(ds10.power, ds11.power),
    connect(ds11.power, ds20.power),
    connect(ds20.power, ds21.power),
    connect(ds21.power, ds30.power),
    connect(ds30.power, ds31.power),
    connect(ds31.power, ds40.power),
    connect(ds40.power, ds41.power),
)
# Impeller assembly
push!(
    cons,
    connect(us40.power, jus.power),
    connect(jus.power, imp.power),
    connect(imp.power, jds.power),
    connect(jds.power, olt.power),
    connect(olt.power, ds10.power),
)
# push!(connect(us30.power, jus.power), connect(jus.power, ds10.power))

# Axis - Impeller coupling
push!(cons, connect(ax.power, agy.pin), connect(agy.pout, imp.power))
push!(cons, agy.Q ~ agy.pin.f, agy.ω ~ agy.pout.f)

# ------------------------------------------------------------------------------
# System

@named psys = ODESystem(cons, t)
mdl = compose(psys, imp, jds, jus, olt, ax, agy, plist...)

generate_graph(mdl)

emdl = expand_connections(mdl)
equations(emdl)
sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)
equations(sys)

# ------------------------------------------------------------------------------
# Human readable variables
@parameters f, fₛ, fᵥ, fₚₒ, fₚᵢ, Iᵢ, Pₜ, Iₚₒ, Iₚᵢ, Iₐ, cₐ, kᵥ, kᵤ, γ_a, γ_b, Qₛ
@variables Qᵢ(t), Qᵤ(t), Qᵥ(t), ω(t), Vᵤ(t), Vᵤ(t), Vᵥ(t)

human = Dict(
    agy.γa => γ_a,
    agy.γb => γ_b,
    imp.m.power.f => Qᵢ,
    ax.m.power.f => ω,
    imp.m.I => Iᵢ,
    agy.ρ => ρ,
    ax.m.I => Iₐ,
    ax.d.R => cₐ * ω,
    ax.T.τ => τ,
)

# Set the friction paramters
human[olt.d.R] = fₛ
human[olt.d.R] = 0
human[imp.d.R] = f

sys = renamevars(sys, human)
equations(sys)
# latexify(sys)

# ------------------------------------------------------------------------------
# Check system

states(sys)
parameters(sys)

# ------------------------------------------------------------------------------
# Simulate

v0 = zeros(length(states(sys)))
rpo["τ"] = 20
vals = Dict(
    Iᵢ => rpo["Iᵢ"],
    γ_b => rpo["γ_b"],
    ρ => rpo["ρ"],
    f => rpo["f"],
    γ_a => rpo["γ_a"],
    cₐ => rpo["cₐ"],
    τ => rpo["τ"],
    Iₐ => rpo["Iₐ"],
    us11.P.p => 0.0,
    ds41.Po.po => 0.0,
)
prob = ODEProblem(sys, v0, (0.0, 30), vals)

sol = solve(prob, reltol = 1e-14, abstol = 1e-14)

# ------------------------------------------------------------------------------
# Plotting

function gen_label(l)
    if length(l) >= 6
        v = ""
        if contains(l, "f")
            v *= "Q"
        else
            v *= "V"
        end

        if contains(l, "us")
            v *= "ᵤ"
        else
            v *= "ᵥ"
        end
        v *= string(l[3]) * "(t)"
        return v
    else
        return l
    end
end

function plotsol(sol, t, idx)
    nsol = sol(t)

    labels = [string(p) for p in nsol.syms[idx]]
    labels = reshape(labels, 1, length(labels))
    n_plots = length(labels)

    plot(nsol.t, nsol[idx, :]', layout = (n_plots, 1), labels = labels, xlabel = "Time (s)")
end

plot(sol.t, sol[us11.plm.power.f])

# plotsol(sol, 0:0.001:10, 3:14)
# plotsol(sol, 0:0.000001:0.0002, 1:2)

# using Measures

# tt = 0:0.00001:10
# nsol = sol(tt)

# idx = 1:14
# labels = [gen_label(string(p)) for p in nsol.syms[idx]]
# labels = reshape(labels, 1, length(labels))
# n_plots = length(labels)

# plot(nsol.t, nsol[idx, :]', layout=(n_plots, 1), titles=labels, size=(800, 140 * n_plots), left_margin=10mm, top_margin=0mm, legend=false, titlefontsize=8)

# nsol.syms

# # ------------------------------------------------------------------------------
# # Filtering

# using PyCall

# pyimport_conda("scipy.optimize", "scipy")
# signal = pyimport("scipy.signal")

# # Aux var for sol with lower discretization
# Δt = 0.001
# tl = 0:Δt:10
# lsol = sol(tl)

# # Filter using bytterworth
# sos = signal.butter(2, 1.5, fs=1/Δt)

# states(sys)
# x = lsol[us11.plm.power.f]
# xf = signal.filtfilt(sos[1], sos[2], x)

# plot(tl, xf, xlabel="Time [s]")
# plot!(tl, xf, xlabel="Time [s]")
# plot!(tl, x)
