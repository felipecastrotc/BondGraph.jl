# # [Duffing Equation](@id msd-duffing)
#
# ## Introduction
#
# The Duffing equation is a nonlinear second-order ordinary differential equation used to model physical systems with nonlinear behaviour, such as electrical circuits, mechanical structures and biological systems. In the mechanical domain it extends the [harmonic oscillator](@ref msd-1dof) with a cubic spring. Despite its simple form, the Duffing equation can show chaotic behaviour. The system is illustrated by the free-body diagram below:
#
# ```@raw html
# <div class="docs-light-only" style="text-align: center;">
# ```
# ![Free-body diagram](../assets/tikz/fbd_duffing.svg)
# ```@raw html
# </div><div class="docs-dark-only" style="text-align: center;">
# ```
# ![Free-body diagram](../assets/tikz/fbd_duffing-dark.svg)
# ```@raw html
# </div>
# ```
#
# ### Duffing oscillator model
#
# We assume a rigid mass $m_{\text{duff}}$ on a linear spring of stiffness $k_{\text{lin}}$, a cubic spring of stiffness $k_{\text{nl}}$ and a viscous damper with coefficient $c_{\text{duff}}$, all massless. The position $x$ is measured from the static equilibrium, positive in the direction of the force $F(t)$.
#
# On the mass act the external force $F(t)$, the linear spring force $-k_{\text{lin}}\,x$, the cubic spring force $-k_{\text{nl}}\,x^3$ and the damper force $-c_{\text{duff}}\,\dot{x}$. Newton's 2nd law gives the governing equation:
#
# ```math
# m_{\text{duff}}\,\ddot{x} + c_{\text{duff}}\,\dot{x} + k_{\text{lin}}\,x + k_{\text{nl}}\,x^3 = F(t).
# ```
#
# With a unit mass and a harmonic force $F(t) = \gamma\cos(\omega t)$ it takes the usual form $\ddot{x} + \delta\,\dot{x} + \alpha\,x + \beta\,x^3 = \gamma\cos(\omega t)$, with $\delta = c_{\text{duff}}$, $\alpha = k_{\text{lin}}$ and $\beta = k_{\text{nl}}$.
#
# See also Kovacic & Brennan (eds.), *The Duffing Equation: Nonlinear Oscillators and their Behaviour*, Wiley 2011, and [Duffing equation](https://en.wikipedia.org/wiki/Duffing_equation).
#
# ### Bond Graph Model
#
# The mass, both springs and the damper move with the same velocity $\dot{x}$, so their bonds meet at a 1-junction, together with the force source:
#
# ```@raw html
# <div class="docs-light-only" style="text-align: center;">
# ```
# ![Bond graph](../assets/tikz/bg_duffing.svg)
# ```@raw html
# </div><div class="docs-dark-only" style="text-align: center;">
# ```
# ![Bond graph](../assets/tikz/bg_duffing-dark.svg)
# ```@raw html
# </div>
# ```
#
# As in the [harmonic oscillator](@ref msd-1dof), the inertia `I` stands for the mass $m_{\text{duff}}$, the compliance `C:1/k_lin` for the linear spring and the resistance `R` for the damper $c_{\text{duff}}$ (Karnopp, Margolis & Rosenberg 2012, §3.1, pp. 38–47). The effort source `Se:F(t)` is the external force. Karnopp's compliance is nonlinear in general, $e = \Phi_C^{-1}(q)$ with $q = \int f\,dt$ (§3.1, pp. 40–42); the cubic spring is such a compliance, `C:1/k_nl`, with $e = k_{\text{nl}}\,q^3$.
#
# ## Preliminary Steps
#
# In this example we first define the cubic compliance as a custom element, then build and simulate the bond graph model, unforced and forced. Finally, we simulate the equation derived above as a plain ModelingToolkit model and compare the results.
#
# We import the bond-graph toolkit to model the system and `OrdinaryDiffEqTsit5` to solve the resulting ODE. We also import ModelingToolkit's independent variable `t` and derivative `D`, for the custom element, the forcing term and the reference model.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using Plots

# ## System Element Definitions
#
# ### Custom compliance
#
# The library's `Compliance` is linear, so the cubic spring is a custom element. `@oneport` adds the element's ports `power` and `pout`, and the element's law is written in its bond variables `e` and `f`, as for the built-in elements: the displacement $q$ integrates the flow, and the effort is $q^3/C$.
@oneport function CubicCompliance(; name, C = 1.0, q = 0.0)
    @parameters C = C           # (m³/N) Cubic compliance, C = 1/k_nl
    @variables q(t) = q         # (m) Displacement, with its initial value
    [D(q) ~ f, e ~ q^3 / C]
end

# ### Elements
#
# The Duffing parameters, in the usual notation:
α = 1.0        # (N/m) Linear stiffness k_lin
β = 5.0        # (N/m³) Cubic stiffness k_nl
δ = 0.02       # (N⋅s/m) Damping coefficient c_duff
γ = 8.0        # (N) Force amplitude
ω = 0.5        # (rad/s) Force frequency

# The mass, the linear and cubic springs and the damper. `Mass`, `Spring` and `Damper` are the mechanical shorthands for `Inertia`, `Compliance` and `Resistance`:
@named mass = Mass(m = 1.0)                     # (kg) Mass m_duff
@named spring = Spring(k = α)                   # (N/m) Linear spring
@named spring_nl = CubicCompliance(C = 1 / β)   # (m³/N) Cubic spring
@named damper = Damper(c = δ)                   # (N⋅s/m) Damper

# The unforced case has a zero force, $F(t) = 0$:
@named force = Se(0.0)                          # (N) External force

# ## Model Construction
#
# A single 1-junction holds the whole model. Following Karnopp's usual orientation (§3.1, pp. 39, 42, 49), the source's bond points into the junction (`force => :pout`, positive) and the bonds of the passive elements point out of it, towards the elements (a bare member `x`). The junction states that all five bonds share the velocity and that the force balance $F - \dot{p} - c_{\text{duff}}\,\dot{x} - k_{\text{lin}}\,q - k_{\text{nl}}\,q_{\text{nl}}^3 = 0$ holds.
@named model = Junction1(force => :pout, mass, damper, spring, spring_nl)

# ## Analysing the model
#
# `mtkcompile` reduces the model to its states: the mass's momentum and the displacements of both springs.
sys = mtkcompile(model)

# Print system states
unknowns(sys)

# Print system parameters
parameters(sys)

# Print the system's equations in terms of the states
full_equations(sys)

# `bg_equations` gives the same equations in terms of the states. We rename the momentum $p$, the springs' displacements $x$ and $x_{\text{nl}}$, and the parameters to the symbols of the hand derivation:
@variables x(t) x_nl(t) p(t)
@parameters m_duff c_duff k_lin k_nl
eqs = bg_equations(model; rename = Dict(
    model.mass.p => p, model.spring.q => x, model.spring_nl.q => x_nl,
    model.mass.I => m_duff, model.damper.R => c_duff,
    model.spring.C => 1 / k_lin, model.spring_nl.C => 1 / k_nl))
latexify(eqs)

# The two displacement equations say both springs' displacements change with the mass's velocity, $\dot{x} = \dot{x}_{\text{nl}} = p/m_{\text{duff}}$, so $x_{\text{nl}} = x$ as long as both start equal; the bond graph keeps them as two states because each compliance stores its own displacement. The momentum equation is the 1-junction's force balance, here with $F = 0$. Substituting $p = m_{\text{duff}}\,\dot{x}$ and $x_{\text{nl}} = x$ gives the Duffing equation derived by hand.
#
# ## System Simulations
#
# ### Bond graph model
#
# Define the simulation time
tspan = (0.0, 10.0)

# The unforced oscillator starts from rest with a unit displacement, $x(0) = 1$, given to both springs.
prob = ODEProblem(sys, [sys.spring.q => 1.0, sys.spring_nl.q => 1.0], tspan)

# Use `solve` to simulate the system.
sol = solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-10)
plot(sol, idxs = [model.spring.q, model.mass.f], labels = ["x" "ẋ"],
     xlabel = "Time (s)", ylabel = "Displacement (m) / velocity (m/s)")

# For the forced case, the force source becomes $F(t) = \gamma\cos(\omega t)$, built from the same elements. The oscillator starts at rest at $x(0) = 0$, the springs' default.
@named force_cos = Se(γ * cos(ω * t))           # (N) External force
@named forced = Junction1(force_cos => :pout, mass, damper, spring, spring_nl)
forced_sol = solve(ODEProblem(mtkcompile(forced), [], tspan), Tsit5(),
                   reltol = 1e-10, abstol = 1e-10)
plot(forced_sol, idxs = [forced.spring.q, forced.mass.f], labels = ["x" "ẋ"],
     xlabel = "Time (s)", ylabel = "Displacement (m) / velocity (m/s)")

# ### Duffing oscillator model
#
# The reference model is the Duffing equation derived by hand, written as two first-order equations in the position $x$ and the velocity $v = \dot{x}$, with a unit mass:
@variables v(t)
duffing_eqs(F) = [D(x) ~ v, D(v) ~ F - δ * v - α * x - β * x^3]

# One `System` per case: the unforced one starts at $x(0) = 1$, the forced one at rest at the origin.
@named ref = System(duffing_eqs(0.0), t)
ref_sol = solve(ODEProblem(mtkcompile(ref), [x => 1.0, v => 0.0], tspan), Tsit5(),
                reltol = 1e-10, abstol = 1e-10)

@named ref_forced = System(duffing_eqs(γ * cos(ω * t)), t)
ref_forced_sol = solve(ODEProblem(mtkcompile(ref_forced), [x => 0.0, v => 0.0], tspan),
                       Tsit5(), reltol = 1e-10, abstol = 1e-10);

# ### Compare models
#
# We compare the bond graph's linear spring displacement with the position $x$ of the reference model.
#
# Unforced case
plot(sol, idxs = model.spring.q, label = "Bond graph", ylabel = "Displacement (m)")
plot!(ref_sol, idxs = x, label = "Model", xlabel = "Time (s)", linestyle = :dash)

# Forced case
plot(forced_sol, idxs = forced.spring.q, label = "Bond graph", ylabel = "Displacement (m)")
plot!(ref_forced_sol, idxs = x, label = "Model", xlabel = "Time (s)", linestyle = :dash)

# The largest difference between the two models on a common time grid, over both springs' displacements and the velocity, for both cases:
tgrid = range(tspan..., length = 1001)
maxdiff(bg, m, ref) = maximum(maximum(abs.(bg(tgrid; idxs = a) .- ref(tgrid; idxs = b)))
                              for (a, b) in [m.spring.q => x, m.spring_nl.q => x,
                                             m.mass.f => v])
maxdiff(sol, model, ref_sol), maxdiff(forced_sol, forced, ref_forced_sol)

# The bond graph model, with its custom cubic compliance, reproduces the hand-derived Duffing dynamics.
