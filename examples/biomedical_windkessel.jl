# # [Arterial Windkessel](@id biomedical-windkessel)
#
# ## Introduction
#
# The heart ejects blood into the aorta only during systole, yet blood flows through the tissues all the time. The large elastic arteries smooth the pulsatile ejection: they stretch and store blood while the heart ejects, and recoil and push it on while the aortic valve is closed. Weber compared this with the air chamber of a fire engine, the *Windkessel*, which turns the pump's strokes into a steady jet, and Frank turned the comparison into a quantitative model (Westerhof, Lankhaar & Westerhof 2009, §1, Fig. 1). The Windkessel is a lumped model: it describes the whole arterial tree by the relation between pressure and flow at its entrance, the aortic root, with a few parameters that have a physiological meaning.
#
# ### Windkessel model
#
# Frank's two-element Windkessel has two parameters (§1):
#
# - the **peripheral resistance** $R$ of the small arteries and arterioles, where most of the resistance to flow is. With venous pressure neglected, it is mean aortic pressure over cardiac output, $R \approx P_{ao,mean}/\mathrm{CO}$;
# - the **total arterial compliance** $C = \Delta V/\Delta P$, the volume the arteries store per unit of pressure, mostly in the large (conduit) arteries.
#
# It predicts that in diastole, with the valve closed and no inflow, aortic pressure decays exponentially, $P(t) = P_{es}\,e^{-t/RC}$ from the end-systolic pressure $P_{es}$ (§2, Fig. 3). It describes systole poorly, though, and its input impedance falls to zero at high frequencies, where the measured one levels off at the **characteristic impedance** $Z_c$ of the proximal aorta. The three-element Windkessel adds $Z_c$ in series at the entrance, before the parallel $C$ and $R$ (§2, Figs. 4–5). $Z_c$ has the dimension of a resistance and is usually drawn as one, although it is not a resistance but a property of wave travel, $Z_c = \mathrm{PWV}\,\rho/A$ (pulse wave velocity, blood density, aortic area); in the systemic circulation it is about 5–7% of $R$ (§2, Fig. 5).
#
# We model the three-element Windkessel as Westerhof et al. draw it (Fig. 5): the aortic inflow $Q(t)$ passes through $Z_c$ into the compliance $C$, which drains through the peripheral resistance $R$ into the veins, whose pressure we take as zero. Let $P$ be the aortic pressure at the entrance and $P_c$ the pressure across the compliance, after $Z_c$. The pressure drop across $Z_c$ is $Z_c Q$, and the inflow splits into the flow that fills the compliance, $C\,\dot{P}_c$, and the peripheral flow $P_c/R$. So
#
# ```math
# C\,\frac{dP_c}{dt} = Q(t) - \frac{P_c}{R}, \qquad P = P_c + Z_c\,Q(t).
# ```
#
# With $Q = 0$ (diastole) $P = P_c$ decays with the time constant $RC$, as for the two-element model. At steady state the mean pressure is $(R + Z_c)$ times the mean flow (§2, Table 1).
#
# The paper doesn't prescribe an inflow; in its examples the inflow is a measured aortic flow (§2, Fig. 6). We assume (our choice, not from the paper) the simplest standard pulsatile waveform: a half-sine during systole and no flow in diastole, with heart period $T$, systole length $T_s$ and stroke volume $SV$,
#
# ```math
# Q(t) = \begin{cases} Q_0 \sin(\pi\,\tau/T_s), & \tau < T_s, \\ 0, & T_s \le \tau < T, \end{cases} \qquad \tau = t \bmod T, \qquad Q_0 = \frac{\pi\,SV}{2\,T_s},
# ```
#
# where $Q_0$ makes the area under one half-sine equal to $SV$.
#
# See also N. Westerhof, J.-W. Lankhaar and B. E. Westerhof, [The arterial Windkessel](https://doi.org/10.1007/s11517-008-0359-2), *Med Biol Eng Comput* 47:131–141, 2009 (open access), Fig. 5 for the hydraulic and electrical forms of the two-, three- and four-element Windkessels.
#
# ### Bond Graph Model
#
# In the hydraulic domain pressure is the effort and volume flow the flow, so their product is power. The heart is a flow source, $Z_c$ and $R$ are resistances and the arterial compliance is a compliance, whose displacement is the stored volume $V = \int Q\,dt$. $Z_c$ carries the whole inflow, so it shares a 1-junction (common flow) with the source; $C$ and $R$ share the pressure $P_c$, so they sit on a 0-junction (common effort). The bond graph is drawn from the model below with `bgplot`, once it is built.
#
# | Element | Stands for | Parameter |
# |---|---|---|
# | `Sf:Q` | the aortic inflow from the heart | $Q(t)$ |
# | `R:Zc` | the characteristic impedance of the proximal aorta | $Z_c$ |
# | `C:C` | the total arterial compliance | $C$ |
# | `R:Rp` | the peripheral resistance | $R$ |
#
# ## Preliminary Steps
#
# In this example we first build and simulate the bond graph model over several heart beats, and check that the pressure settles into a periodic range of plausible values. Then we simulate the hand-derived equation as a plain ModelingToolkit model and compare the results.
#
# We import the bond-graph toolkit to model the system and `OrdinaryDiffEqTsit5` to solve the resulting ODE. We also import ModelingToolkit's independent variable `t` and derivative `D`, for the inflow and the reference model.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using Plots

# ## System Element Definitions
#
# We work in the units of clinical hemodynamics: pressure in mmHg, volume in mL, time in s, so resistances are in mmHg⋅s/mL and the compliance in mL/mmHg. Westerhof et al. (2009) give no numerical parameter values, only that $Z_c$ is about 5–7% of $R$ (§2). The values below are our assumption, chosen to give an adult at rest a plausible mean pressure and cardiac output.
#
# ### Aortic inflow
#
# A heart rate of 75 beats/min, a systole of 0.3 s and a stroke volume of 70 mL, that is a cardiac output of 5.25 L/min:
T = 0.8                         # (s) Heart period
Ts = 0.3                        # (s) Duration of systole (ejection)
SV = 70.0                       # (mL) Stroke volume
Q0 = π * SV / (2Ts)             # (mL/s) Peak inflow

# The half-sine waveform, as a function of time. Called with the symbolic `t` it gives the flow source its law; we use the same function in the reference model.
inflow(t) = ifelse(mod(t, T) < Ts, Q0 * sin(π * mod(t, T) / Ts), 0.0)

# ### Elements
#
# The peripheral resistance gives a mean pressure of about $R\,SV/T \approx 88$ mmHg; $Z_c$ is 6% of it, in the middle of the paper's range; the compliance gives a diastolic decay time $RC = 1.5$ s, about twice the heart period:
@named Q = Sf(inflow(t))                    # (mL/s) Aortic inflow
@named Zc = Resistance(R = 0.06)            # (mmHg⋅s/mL) Characteristic impedance
@named Rp = Resistance(R = 1.0)             # (mmHg⋅s/mL) Peripheral resistance
@named C = Compliance(C = 1.5, q = 120.0)   # (mL/mmHg, mL) Total arterial compliance

# The compliance's state is the volume it stores above its unstressed volume. Its initial value, 120 mL, starts the pressure $P_c = V/C$ at 80 mmHg, a typical diastolic pressure.
#
# ## Model Construction
#
# The inflow, $Z_c$ and the bond to the rest of the tree share one flow, $Q$: a 1-junction. The source's bond points into the junction, away from the source, and $Z_c$'s bond points out of the junction into it (`[-1, Zc]`), Karnopp, Margolis & Rosenberg's usual orientations (2012, §3.1, pp. 39 and 49). The bond to the compliance and the periphery leaves the junction (`[-1, :wk]`), following the blood downstream. The effort sum is therefore $P - Z_c Q - P_c = 0$: the source's effort, the pressure at the entrance, is $P = Z_c Q + P_c$.
@named j1 = Junction1(Q, [-1, Zc], [-1, :wk])

# $C$ and $R$ share the pressure $P_c$: a 0-junction. The bond from the 1-junction enters it (`:wk`), and the bonds to $C$ and $R$ point out of it into the elements (`[-1, C]`, `[-1, Rp]`), again the usual orientations. The flow balance is $Q - \dot{V} - P_c/R = 0$, the volume balance of the arteries. The veins are the pressure reference, so $R$'s bond ends there with no element.
@named j0 = Junction0(:wk, [-1, C], [-1, Rp])

# The two free ports `wk` are the two ends of one bond, joined with `connect`:
@named windkessel = System([connect(j1.wk, j0.wk)], t; systems = [j1, j0])

# The library draws the bond graph from the model, in the order of Westerhof et al.'s Fig. 5: source, $Z_c$, then $C$ and $R$:
bgplot(windkessel)

# ## Analysing the model
#
# `mtkcompile` reduces the model to its one state: the volume stored in the compliance.
sys = mtkcompile(windkessel)

# Print system states
unknowns(sys)

# Print system parameters
parameters(sys)

# Print the system's equations in terms of the states
full_equations(sys)

# ### Generated equations against the hand derivation
#
# `bg_equations` returns the same equations with the names we choose. We rename the stored volume to $V$ and the parameters to the symbols of the introduction (inside a `let`, so the name `C` keeps pointing at the element):
let
    @variables V(t)
    @parameters C R
    bg_equations(windkessel; rename = Dict(j0.C.q => V, j0.C.C => C, j0.Rp.R => R))
end

# The inflow appears as its expression in $t$. With $V = C\,P_c$, the equation reads $C\,\dot{P}_c = Q(t) - P_c/R$, the hand-derived equation. $Z_c$ doesn't appear, and renaming it would be an error: it carries the source's flow, which is known, so it only sets the source's effort, the entrance pressure $P = P_c + Z_c\,Q$, an observed variable that `bg_equations` (state equations only) doesn't list. We check it in the simulation below.
#
# ## System Simulations
#
# ### Bond graph model
#
# Define the simulation time: ten heart beats.
tspan = (0.0, 10T)

# Generate an `ODEProblem`, with the initial volume given to the compliance:
prob = ODEProblem(sys, [], tspan)

# Use `solve` to simulate the system. Tight tolerances let us compare it with the reference model at the end. The inflow's slope jumps at the start and end of each systole; `tstops` puts a solver step on each of these instants.
kinks = sort([(0:9) * T; (0:9) * T .+ Ts])
sol = solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-10, tstops = kinks)

# The entrance (aortic) pressure is the source's effort, and the pressure after $Z_c$ the compliance's:
plot(sol, idxs = [j1.Q.power.e, j0.C.power.e], labels = ["P (aortic)" "P_c"],
     xlabel = "Time (s)", ylabel = "Pressure (mmHg)")

# The inflow, the source's flow, against the peripheral outflow, $R$'s flow. The compliance turns the pulsatile inflow into a nearly steady outflow:
plot(sol, idxs = [j1.Q.power.f, j0.Rp.power.f], labels = ["Q (aortic inflow)" "P_c/R (peripheral)"],
     xlabel = "Time (s)", ylabel = "Flow (mL/s)")

# ### Periodic pressure range
#
# After a few beats the pressure repeats from beat to beat. The diastolic (minimum) and systolic (maximum) aortic pressure of each beat, on a fine grid:
beats = [range(k * T, (k + 1) * T, length = 801) for k in 0:9]
P_range = [extrema(sol(b; idxs = j1.Q.power.e).u) for b in beats]

# The last beat's diastolic and systolic pressures, about 73/117 mmHg, and their change from the beat before, a few hundredths of a mmHg: the start-up transient decays with the time constant $RC = 1.5$ s, by a factor $e^{-T/RC} \approx 0.6$ per beat, and after ten beats is almost gone:
P_range[end], P_range[end] .- P_range[end-1]

# The mean aortic pressure over the last beat and the predicted $(R + Z_c)\,SV/T$:
let b = beats[end]
    sum(sol(b[1:end-1]; idxs = j1.Q.power.e).u) / (length(b) - 1), (1.0 + 0.06) * SV / T
end

# ### Windkessel model
#
# The reference model is the hand-derived equation, in the pressure $P_c$, with the same parameter values and initial pressure:
Rᵥ, Zcᵥ, Cᵥ = 1.0, 0.06, 1.5     # (mmHg⋅s/mL, mmHg⋅s/mL, mL/mmHg)
@variables Pc(t) = 80.0          # (mmHg) Pressure across the compliance
@variables P(t)                  # (mmHg) Aortic pressure at the entrance
@named ref = System([Cᵥ * D(Pc) ~ inflow(t) - Pc / Rᵥ, P ~ Pc + Zcᵥ * inflow(t)], t)

# We compile and solve it as the bond graph model, with the same time span, tolerances and stops.
ref_sol = solve(ODEProblem(mtkcompile(ref), [], tspan), Tsit5(), reltol = 1e-10, abstol = 1e-10,
                tstops = kinks);

# ### Compare models
#
# The aortic pressure of the bond graph (solid) against the Windkessel equation (dashed):
plot(sol, idxs = j1.Q.power.e, label = "Bond graph", xlabel = "Time (s)", ylabel = "Aortic pressure (mmHg)")
plot!(ref_sol, idxs = P, label = "Windkessel equation", linestyle = :dash)

# On a common time grid, the largest difference between the two models, in the aortic pressure $P$ and the pressure after $Z_c$, $P_c$:
tgrid = range(tspan..., length = 1001)
(maximum(abs, sol(tgrid; idxs = j1.Q.power.e).u - ref_sol(tgrid; idxs = P).u),
 maximum(abs, sol(tgrid; idxs = j0.C.power.e).u - ref_sol(tgrid; idxs = Pc).u))

# The difference is about $10^{-9}$ of the pressure, the solver's tolerance. The bond graph model reproduces the three-element Windkessel, and its aortic pressure settles into a periodic beat with plausible systolic and diastolic values.
