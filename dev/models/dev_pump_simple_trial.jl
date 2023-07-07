using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations
using ModelingToolkit

# ==========================================================
# Used alias
@variables ω(t), Q(t)

# ==========================================================
# New test with multiple outputs for pump leakage

# ----------------------------------------------------------
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
@named plm = Mass(m=rpo["Iₚᵢ"])
@named pld = Damper(c=rpo["fₚᵢ"])     # General damping element
@named pls = Spring(k=-1 / rpo["kᵤ"])     # General compliance element

# ----------------------------------------------------------
# Impeller Junctions

@named lek = Junction1([-1, d])     # Leakage setting

@named jus = Junction0()     # Coupling for leakage

@named imp = Junction1([-1, m], [-1, d])     # Impeller

@named olt = Junction1([-1, d])     # Shock losses
@named jds = Junction0()            # Coupling for leakage

# ----------------------------------------------------------
# Axis Junctions
@named ax = Junction1(T, [-1, m], [-1, d])      # Axis parameters

# ----------------------------------------------------------
# Pipe Junctions

@named us11 = Junction1([-1, plm], [-1, pld], P)     # Inlet pipe
@named us10 = Junction0([-1, pls])
# @named us10 = Junction0()
@named us21 = Junction1([-1, plm], [-1, pld])     # Inlet pipe
@named us20 = Junction0([-1, pls])
# @named us20 = Junction0()
@named us31 = Junction1([-1, plm], [-1, pld])     # Inlet pipe
@named us30 = Junction0([-1, pls])
# @named us30 = Junction0()
plist = [us11, us10, us21, us20, us31, us30];

@named ds10 = Junction0([-1, pls])
# @named ds10 = Junction0()
@named ds11 = Junction1([-1, plm], [-1, pld])     # Outlet pipe
@named ds20 = Junction0([-1, pls])
# @named ds20 = Junction0()
@named ds21 = Junction1([-1, plm], [-1, pld])     # Outlet pipe
@named ds30 = Junction0([-1, pls])
# @named ds30 = Junction0()
@named ds31 = Junction1([-1, plm], [-1, pld], Po)     # Outlet pipe
push!(plist, ds11, ds10, ds21, ds20, ds31, ds30);

# Tank junction
@named tnkin = Junction0(P)
@named tnkout = Junction0(Po)

# ----------------------------------------------------------
# MGY
@parameters γa, γb

gyp = ρ * (γa * ω - γb * Q)
@named agy = mGY(g=gyp)

# ----------------------------------------------------------
# Build connections

# Impeller
# cons = [connect(ilt.power, jus.power), connect(lek.power, jus.power), connect(jus.power, imp.power), connect(imp.power, jds.power), connect(jds.power, lek.power), connect(jds.power, olt.power)]# With impeller leakage
# cons = [connect(ilt.power, jus.power), connect(jus.power, imp.power), connect(imp.power, jds.power), connect(jds.power, olt.power)] # Without impeller leakage
# cons = [connect(tnk.power, ilt.power), connect(ilt.power, jus.power), connect(jus.power, imp.power), connect(imp.power, jds.power), connect(jds.power, olt.power)] # Without impeller leakage
# cons = [connect(tnk.power, ilt.power), connect(ilt.power, ilt0.power), connect(ilt0.power, jus.power), connect(jus.power, imp.power), connect(imp.power, jds.power), connect(jds.power, olt.power)] # Without impeller leakage
# cons = [connect(tnk.power, ilt.power), connect(ilt.power, jus.power), connect(jus.power, imp.power), connect(imp.power, jds.power), connect(jds.power, olt.power)] # Without impeller leakage

# Pipeline assembly
cons = [connect(us11.power, us10.power), connect(us10.power, us21.power), connect(us21.power, us20.power), connect(us20.power, us31.power), connect(us31.power, us30.power)]
push!(cons, connect(ds10.power, ds11.power), connect(ds11.power, ds20.power), connect(ds20.power, ds21.power), connect(ds21.power, ds30.power), connect(ds30.power, ds31.power))
# Impeller assembly
push!(cons, connect(us30.power, jus.power), connect(jus.power, imp.power), connect(imp.power, jds.power), connect(jds.power, olt.power), connect(olt.power, ds10.power))
# push!(connect(us30.power, jus.power), connect(jus.power, ds10.power))

# Axis - Impeller coupling
push!(cons, connect(ax.power, agy.pin), connect(agy.pout, imp.power))
push!(cons, agy.Q ~ agy.pin.f, agy.ω ~ agy.pout.f)

# ----------------------------------------------------------
# System

@named psys = ODESystem(cons, t)
# mdl = compose(psys, imp, jds, jus, ilt, olt, ax, agy)
# mdl = compose(psys, lek, imp, jds, jus, ilt, olt, ax, agy)
# mdl = compose(psys, lek, imp, jds, jus, ilt, olt, ax, agy, plds)
# mdl = compose(psys, lek, imp, jds, jus, ilt, olt, ax, agy, plds, plus)
# mdl = compose(psys, imp, jds, jus, ilt, olt, ax, agy, plds)
# mdl = compose(psys, lek, imp, jds, jus, ilt, olt, plds, plus)
# mdl = compose(psys, imp, jds, jus, ilt, olt, plds, plus)
# mdl = compose(psys, imp, jds, jus, ilt, olt, plds, tnk, ilt0, plds0, ax, agy)
# mdl = compose(psys, imp, jds, jus, ilt, olt, plds, tnk, ilt0, plds0, tnk2)
# mdl = compose(psys, imp, jds, jus, ilt, olt, plds, tnk, tnk2, ax, agy)
mdl = compose(psys, imp, jds, jus, olt, ax, agy, plist...)
# mdl = compose(psys, plist..., jus)


generate_graph(mdl)

emdl = expand_connections(mdl)
equations(emdl)

amdl = ModelingToolkit.alias_elimination(emdl)

equations(amdl)
sys = ModelingToolkit.structural_simplify(emdl)
@named sys = reducedobs(sys)

equations(sys)

@parameters f, fₛ, fᵥ, fₚₒ, fₚᵢ, Iᵢ, Pₜ, Iₚₒ, Iₚᵢ, Iₐ, cₐ, kᵥ, kᵤ, γ_a, γ_b, Qₛ
@variables Qᵢ(t), Qᵤ(t), Qᵥ(t), ω(t), Vᵤ(t), Vᵤ(t), Vᵥ(t)

# olt.d.R => fₛ / Qᵢ * (1 - Qᵢ / Qₛ)^2
# human = Dict(agy.γa => γ_a, agy.γb => γ_b, imp.m.power.f => Qᵢ, ax.m.power.f => ω, imp.m.I => Iᵢ, tnk.P.p => Pₜ, plds.m2.power.f => Qᵥ, plds.m2.I => Iₚₒ, ilt.m.I => Iₚᵢ, agy.ρ => ρ, ax.m.I => Iₐ, ax.d.R => cₐ * ω, ilt.m.power.f => Qᵤ, ax.T.τ => τ, ilt0.s.q => Vᵤ, plds0.s.q => Vᵥ, plds0.s.C => kᵥ, ilt0.s.C => kᵤ, tnk2.P.p => Pₜ)
human = Dict(agy.γa => γ_a, agy.γb => γ_b, imp.m.power.f => Qᵢ, ax.m.power.f => ω, imp.m.I => Iᵢ, agy.ρ => ρ, ax.m.I => Iₐ, ax.d.R => cₐ * ω, ax.T.τ => τ)

# Set the friction paramters
human[olt.d.R] = fₛ
human[olt.d.R] = 0
human[imp.d.R] = f

sys = renamevars(sys, human)

equations(sys)
# latexify(sys)

states(sys)
parameters(sys)

# param_range = Dict(
#     "Iᵢ" => Iᵢ,
#     "Iₐ" => Iₐ,
#     "Iₚᵢ" => I_l,
#     "Iₚₒ" => I_l,
#     "f" => fₓ,
#     "fₛ" => kₛ,
#     "Qₛ" => Qₛ,
#     "fₚᵢ" => fₓ,
#     "fₚₒ" => fₓ,
#     "fᵥ" => fₓ,
#     "kᵥ" => k,
#     "kᵤ" => k,
#     "cₐ" => cₐ,
#     "γ_a" => γ_a,
#     "γ_b" => γ_b,
#     "Pₜ" => [0, 0],
#     "τ" => τ,
# )

states(sys)

param_range["τ"][2] = 1

v0 = [0.0, 0.0, 0.0, 0.0, V₀[2], V₀[2]]
v0 = zeros(length(states(sys)))
# v0 = [0.0, 0.0, 0.0]
# rpo["τ"] = -0.01
# rpo["γ_b"] = 3.3
# rpo["γ_a"] = 0.02
# rpo["Pₜ"]
rpo["τ"] = 10
#vals = Dict(p => float(param_range[string(p)][2]) for p in parameters(sys))
# vals = Dict(p => float(rpo[string(p)]) for p in parameters(sys))
vals = Dict(Iᵢ => rpo["Iᵢ"], γ_b => rpo["γ_b"], ρ => rpo["ρ"], f => rpo["f"], γ_a => rpo["γ_a"], cₐ => rpo["cₐ"], τ => rpo["τ"], Iₐ => rpo["Iₐ"], us11.P.p => 0.0, ds31.Po.po => 0.0)
prob = ODEProblem(sys, v0, (0.0, 10), vals)
sol = solve(prob, reltol=1e-14, abstol=1e-14)
# sol = solve(prob)

states(sys)
plot(sol.t, sol[Qᵢ])
plot(sol.t, sol[us10.pls.q])
plot(sol)

plot(sol.t, sol[us31.plm.power.f], label="Qu")
plot!(sol.t, sol[ds11.plm.power.f], label="Qv")

plot(sol.t, sol[ω], label="Qu")
plot!(sol.t, sol[Qᵢ], label="Qv")

nsol = sol(0:0.00001:0.0005);
plot(nsol.t, [nsol[ω], nsol[Qᵢ]], layout=(2, 1));

plot(sol, layout=(14, 1))

equations(sys)
states(sys)

# ------------------------------------------------------------------------------

l = @layout [grid(3, 1)]
p1 = plot(...)
p2 = plot(...)
p3 = plot(...)
plot(p1, p2, p3, layout=l)


tv = 0:0.001:10
nnsol = sol(tv)
q = nsol[ds31.plm.power.f]
q = nsol[us31.plm.power.f]

t_fast = 0:0.000001:0.0005
sol_fast = sol(t_fast)

labels = reshape(v, length(v), 1)

[string(p) for p in sol_fast.syms[1:2]][:, :]

plot(t_fast, sol_fast, layout=(2, 1), label=[])

plot(sol_fast.t, sol_fast[1:4, :]', layout=(4, 1), labels = tits)

plot(rand(100, 4), layout=(4), label=["a" "b" "c" "d"],
    title=["1" "2" "3" "4"])

string.(sol_fast.syms[1:2])
length(sol_fast.t)

size(sol_fast)




# ------------------------------------------------------------------------------
# Filtering
using PyCall
pyimport_conda("scipy.optimize", "scipy")
signal = pyimport("scipy.signal")

sos = signal.butter(5, 1.5, fs=1000)
filtered = signal.filtfilt(sos[1], sos[2], q)

plot(tv, filtered)
plot!(tv, q)