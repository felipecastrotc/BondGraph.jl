# # [Coupled Tanks](@id hydraulic-coupled-tanks)
#
# ## Introduction
#
# Two tanks in series are the standard liquid-level process plant: a pump fills the first tank, the first tank drains through a valve into the second, and the second drains through a valve to the atmosphere. The valves are orifices, whose flow grows with the square root of the pressure drop across them, so the model is nonlinear. Borutzky (2011, Ch. 4, §4.4.1, pp. 145–148, Fig. 4.7) uses this plant as an example; a schematic is in his Fig. 4.7.
#
# ### Coupled tanks model
#
# We assume an incompressible liquid of density $\rho$, open tanks with vertical walls and cross-section areas $A_1$, $A_2$, and outlets at the bottom of each tank. The pressure at the bottom of tank $i$, above the atmospheric pressure, is the hydrostatic pressure of its liquid column of level $h_i$ and volume $V_i = A_i h_i$:
#
# ```math
# p_i = \rho\,g\,h_i = \frac{\rho\,g}{A_i}\,V_i .
# ```
#
# A tank is therefore a linear hydraulic capacitance, $V_i = C_i\,p_i$ with $C_i = A_i/(\rho g)$.
#
# Each valve is an orifice with Bernoulli's square-root law, as Borutzky writes it (Ch. 4, §4.2.1, p. 141, Eq. 4.12):
#
# ```math
# Q = c_d\,a\,\mathrm{sign}(\Delta p)\,\sqrt{\frac{2\,|\Delta p|}{\rho}},
# ```
#
# with discharge coefficient $c_d$, orifice area $a$ and pressure drop $\Delta p$ in the direction of $Q$. With $\Delta p = \rho g \Delta h$ it is Torricelli's law, $Q = c_d\,a\,\mathrm{sign}(\Delta h)\sqrt{2 g |\Delta h|}$. The sign-and-absolute-value form gives the flow for either direction of the pressure drop. The plain $\sqrt{\Delta p}$ is undefined (a `DomainError` in Julia) as soon as $\Delta p < 0$: that happens physically when the second tank is fuller than the first and liquid flows back through valve 1, and numerically whenever a solver step overshoots a level that is close to zero.
#
# The volume balance of each tank, with the pump flow $Q_p$, the flow $Q_1$ through valve 1 and the outflow $Q_o$ through valve 2, gives the governing equations in the levels:
#
# ```math
# \begin{aligned}
# A_1\,\dot{h}_1 &= Q_p - Q_1, & Q_1 &= c_{d}\,a_1\,\mathrm{sign}(h_1 - h_2)\sqrt{2 g\,|h_1 - h_2|},\\
# A_2\,\dot{h}_2 &= Q_1 - Q_o, & Q_o &= c_{d}\,a_2\,\mathrm{sign}(h_2)\sqrt{2 g\,|h_2|}.
# \end{aligned}
# ```
#
# At steady state both flows equal the pump flow, $Q_1 = Q_o = Q_p$, so the levels are
#
# ```math
# h_2^* = \frac{1}{2g}\left(\frac{Q_p}{c_d\,a_2}\right)^2, \qquad
# h_1^* = h_2^* + \frac{1}{2g}\left(\frac{Q_p}{c_d\,a_1}\right)^2 .
# ```
#
# Borutzky linearises the valves around an operating point to get a linear state-space model in the tank pressures (Eq. 4.28), for his parameter-sensitivity analysis. We keep the nonlinear orifices: that is our modelling choice, not the book's. His analytical redundancy relation for the second tank uses the nonlinear valve law in the same sign-safe form, $k_1\,\mathrm{sign}(p_1 - p_2)\sqrt{|p_1 - p_2|}$ (Ch. 4, §4.6.1, p. 167, Eq. 4.94).
#
# See also W. Borutzky, *Incremental Bond Graphs*, Ch. 4 of W. Borutzky (ed.), *Bond Graph Modelling of Engineering Systems*, Springer 2011, Figs. 4.7–4.8.
#
# ### Bond Graph Model
#
# Pressure is the effort and volume flow the flow, so their product is power. Each tank's pressure is common to the flows in and out of it, so each tank is a C element on a 0-junction (common effort). Each valve carries one flow driven by the difference of the pressures on its two sides, so it is an R element on a 1-junction (common flow, efforts summed). The pump is a flow source, and the atmosphere is an effort source of zero (gauge) pressure. This is the bond graph of Borutzky's Fig. 4.8 (p. 146), without its detectors. The library draws it from the model below with `bgplot`, once it is built.
#
# | Element | Stands for | Parameter |
# |---|---|---|
# | `Sf:Qp` | the pump | $Q_p$ |
# | `C:tank1`, `C:tank2` | the tanks | $C_i = A_i/(\rho g)$ |
# | `R:valve1`, `R:valve2` (drawn `Orifice:valve1`, `Orifice:valve2`) | the valves (orifices) | $c_d$, $a_i$, $\rho$ |
# | `Se:atm` | the atmosphere | $p = 0$ |
#
# ## Preliminary Steps
#
# In this example we first define the orifice as a custom element, then build and simulate the bond graph model and check its steady state against the closed form. Finally, we simulate the level equations derived above as a plain ModelingToolkit model and compare the results.
#
# We import the bond-graph toolkit to model the system and `OrdinaryDiffEqTsit5` to solve the resulting ODE. We also import ModelingToolkit's independent variable `t` and derivative `D`, for the reference model.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using Plots

# ## System Element Definitions
#
# ### Custom orifice
#
# The library's `Resistance` is linear, so the orifice is a custom element. `@oneport` adds the element's power port `power`, and the law is written in its bond variables, Eq. 4.12 with the pressure drop `power.e` and the flow `power.f`. It is in conductance form, the flow explicit in the effort:
@oneport function Orifice(; name, cd = 0.6, a = 1e-4, ρ = 1000.0)
    @parameters cd = cd         # (-) Discharge coefficient
    @parameters a = a           # (m²) Orifice area
    @parameters ρ = ρ           # (kg/m³) Liquid density
    [power.f ~ cd * a * sign(power.e) * sqrt(2 * abs(power.e) / ρ)]
end

# ### Parameters
#
# Borutzky gives no numerical values for this example. Ours are a water plant of laboratory size: two tanks of 0.1 m² cross-section, valves of 5 cm² with a sharp-edged orifice's discharge coefficient, and a pump delivering 1 L/s.
ρ = 1000.0          # (kg/m³) Water density
g = 9.81            # (m/s²) Gravitational acceleration
A1 = 0.1            # (m²) Tank 1 cross-section
A2 = 0.1            # (m²) Tank 2 cross-section
a1 = 5e-4           # (m²) Valve 1 orifice area
a2 = 5e-4           # (m²) Valve 2 orifice area
cd = 0.6            # (-) Discharge coefficient of both valves
Qp = 1e-3           # (m³/s) Pump flow

# Tank 1 starts empty and tank 2 at 1 m, above both steady-state levels, so liquid first flows back from tank 2 into tank 1 and the pressure drop across valve 1 changes sign:
h1₀, h2₀ = 0.0, 1.0  # (m) Initial levels

# ### Elements
#
# The tanks are compliances with $C = A/(\rho g)$ and the stored volume as state, $q = V = A h$:
@named tank1 = Compliance(C = A1 / (ρ * g), q = A1 * h1₀)   # (m³/Pa, m³) Tank 1
@named tank2 = Compliance(C = A2 / (ρ * g), q = A2 * h2₀)   # (m³/Pa, m³) Tank 2

# The valves, the pump and the atmosphere:
@named valve1 = Orifice(cd = cd, a = a1, ρ = ρ)   # Valve between the tanks
@named valve2 = Orifice(cd = cd, a = a2, ρ = ρ)   # Outlet valve
@named pump = Sf(Qp)                               # (m³/s) Pump flow Qp
@named atm = Se(0.0)                               # (Pa) Atmospheric (gauge) pressure

# ## Model Construction
#
# The bonds follow Borutzky's Fig. 4.8: every half-arrow points downstream, from the pump to the atmosphere, so a positive flow on any bond is liquid moving towards the outlet.
#
# Tank 1's pressure $p_1$ is common to the pump, the tank and valve 1: a 0-junction. The pump's bond points into it, away from the source; the tank's bond points out of it into the tank (`[-1, tank1]`), so the tank fills when its flow is positive; the bond to valve 1 leaves it (`[-1, :v1]`). The flow balance is $Q_p - \dot{V}_1 - Q_1 = 0$.
@named j1 = Junction0(pump, [-1, tank1], [-1, :v1])

# Valve 1 carries one flow $Q_1$ from tank 1 to tank 2: a 1-junction. The bond from tank 1 enters it (`:v1`), the bond to the valve points into the valve (`[-1, valve1]`) and the bond to tank 2 leaves it (`[-1, :t2]`). The effort sum is $p_1 - \Delta p_1 - p_2 = 0$, so the valve sees the drop $\Delta p_1 = p_1 - p_2$, and the drop and the flow are both positive from tank 1 to tank 2.
@named jv1 = Junction1(:v1, [-1, valve1], [-1, :t2])

# Tank 2 is a 0-junction like tank 1, filled by valve 1's bond (`:t2`) and drained by the bond to valve 2 (`[-1, :v2]`): $Q_1 - \dot{V}_2 - Q_o = 0$.
@named j2 = Junction0(:t2, [-1, tank2], [-1, :v2])

# Valve 2 carries the outflow $Q_o$ from tank 2 to the atmosphere. The bond to the atmosphere points into the `Se:0` element (`[-1, atm]`), as in the book's figure: the atmosphere receives the outflow. The effort sum $p_2 - \Delta p_2 - 0 = 0$ makes the outlet drop $\Delta p_2 = p_2$. The order of the members only sets the drawing: listed first, the atmosphere is drawn on the power path, at the right end as in the book.
@named jv2 = Junction1(:v2, [-1, atm], [-1, valve2])

# Each pair of free ports with the same name is one bond, joined with `connect`:
@named tanks = System([connect(j1.v1, jv1.v1), connect(jv1.t2, j2.t2), connect(j2.v2, jv2.v2)], t;
                      systems = [j1, jv1, j2, jv2])

# The library draws the bond graph from the model, in the order of Borutzky's Fig. 4.8. It labels each element with the constructor that built it, so the valves, R elements in the book (`R : R1()`), appear as `Orifice`:
bgplot(tanks)

# ## Analysing the model
#
# `mtkcompile` reduces the model to its states: the volumes stored in the two tanks.
sys = mtkcompile(tanks)

# Print system states
unknowns(sys)

# Print system parameters
parameters(sys)

# Print the system's equations in terms of the states
full_equations(sys)

# ### Generated equations against the hand derivation
#
# `bg_equations` returns the same equations with the names we choose. We rename the stored volumes to $V_1$, $V_2$ and the parameters to symbols of the introduction (inside a `let`, so the Julia values `A1`, `cd`, … keep their numbers):
let
    @variables V1(t) V2(t)
    @parameters C1 C2 cd a1 a2 ρ
    bg_equations(tanks; rename = Dict(j1.tank1.q => V1, j2.tank2.q => V2,
                                      j1.tank1.C => C1, j2.tank2.C => C2,
                                      jv1.valve1.cd => cd, jv2.valve2.cd => cd,
                                      jv1.valve1.a => a1, jv2.valve2.a => a2,
                                      jv1.valve1.ρ => ρ, jv2.valve2.ρ => ρ))
end

# The pump flow appears as its value. The tank pressures are $p_i = V_i/C_i$, so the arguments of the square roots are $p_1 - p_2$ for valve 1 and $p_2$ for valve 2, as the two 1-junctions require. With $V_i = A_i h_i$, $C_i = A_i/(\rho g)$ and $p_i = \rho g h_i$, and $\sqrt{2\rho g |\Delta h|/\rho} = \sqrt{2 g |\Delta h|}$, the two equations are the level equations of the introduction.
#
# ## System Simulations
#
# ### Bond graph model
#
# Define the simulation time: long enough for the levels to settle.
tspan = (0.0, 5000.0)

# Generate an `ODEProblem`, with the initial volumes given to the tanks:
prob = ODEProblem(sys, [], tspan)

# Use `solve` to simulate the system. Tight tolerances let us compare it with the reference model at the end.
sol = solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-12)

# The levels are the stored volumes over the tank areas:
plot(sol, idxs = [j1.tank1.q / A1, j2.tank2.q / A2], labels = ["h₁" "h₂"],
     xlabel = "Time (s)", ylabel = "Level (m)")

# The flow through valve 1 is negative at first, liquid flowing back from tank 2 into tank 1, until the levels cross; both valve flows then rise to the pump flow:
plot(sol, idxs = [jv1.valve1.power.f, jv2.valve2.power.f], labels = ["Q₁ (valve 1)" "Qₒ (valve 2)"],
     xlabel = "Time (s)", ylabel = "Flow (m³/s)", xlims = (0, 1000))

# ### Steady state
#
# The closed-form steady-state levels of the introduction:
h2ss = (Qp / (cd * a2))^2 / 2g
h1ss = h2ss + (Qp / (cd * a1))^2 / 2g
(h1ss, h2ss)

# The difference between the simulated levels at the end of the run and the closed form, small once the transient has decayed:
(sol[j1.tank1.q][end] / A1 - h1ss, sol[j2.tank2.q][end] / A2 - h2ss)

# ### Coupled tanks model
#
# The reference model is the level equations derived by hand, with the same parameter values and initial levels. It is written in Torricelli's form, with $g$ and no density, so it does not reuse the bond graph's pressures:
@variables h1(t) = h1₀ h2(t) = h2₀   # (m) Tank levels
@variables Q1(t) Qo(t)               # (m³/s) Valve flows
@named ref = System([A1 * D(h1) ~ Qp - Q1,
                     A2 * D(h2) ~ Q1 - Qo,
                     Q1 ~ cd * a1 * sign(h1 - h2) * sqrt(2g * abs(h1 - h2)),
                     Qo ~ cd * a2 * sign(h2) * sqrt(2g * abs(h2))], t)

# We compile and solve it as the bond graph model, with the same time span and tolerances.
ref_sol = solve(ODEProblem(mtkcompile(ref), [], tspan), Tsit5(), reltol = 1e-10, abstol = 1e-12);

# ### Compare models
#
# The levels of the bond graph (solid) against the level equations (dashed):
plot(sol, idxs = [j1.tank1.q / A1, j2.tank2.q / A2], labels = ["h₁ bond graph" "h₂ bond graph"],
     xlabel = "Time (s)", ylabel = "Level (m)")
plot!(ref_sol, idxs = [h1, h2], labels = ["h₁ level equations" "h₂ level equations"], linestyle = :dash)

# On a common time grid, the largest difference between the two models, in the two levels:
tgrid = range(tspan..., length = 1001)
(maximum(abs, sol(tgrid; idxs = j1.tank1.q / A1).u - ref_sol(tgrid; idxs = h1).u),
 maximum(abs, sol(tgrid; idxs = j2.tank2.q / A2).u - ref_sol(tgrid; idxs = h2).u))

# The bond graph model, with its custom sign-safe orifice, reproduces the hand-derived level dynamics, through the reversal of the flow in valve 1, and settles at the closed-form steady-state levels.
