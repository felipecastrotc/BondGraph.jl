# # [Chemical Reaction Network](@id chemical-reaction-network)
#
# ## Introduction
#
# A chemical reaction network converts species into one another. Its thermodynamics fixes where it can go: every reaction runs down a difference of chemical potential, and at equilibrium all those differences vanish. Gawthrop & Crampin (2014) model such networks as bond graphs, with chemical potential as the effort and molar flow as the flow, so that the models obey the laws of thermodynamics by construction. This example builds their two-reaction chain (Gawthrop & Crampin 2014, Fig. 1c):
#
# ```math
# \mathrm{A} \rightleftharpoons \mathrm{B} \rightleftharpoons \mathrm{C}.
# ```
#
# ### Reaction network model
#
# We assume a dilute, isothermal and isobaric mixture in a fixed volume, so that each species has the chemical potential (Gawthrop & Crampin 2014, Eqs. 2–3)
#
# ```math
# \mu_i = \mu_i^0 + RT \ln\frac{x_i}{V} = RT \ln(K_i\,x_i), \qquad K_i = \frac{1}{V}\,e^{\mu_i^0/RT},
# ```
#
# where $x_i$ is the amount of species $i$ (mol), $\mu_i^0$ its standard chemical potential, $R$ the gas constant and $T$ the temperature. The thermodynamic constant $K_i$ (1/mol) gathers the standard potential and the volume.
#
# Each reaction's molar flow follows the Marcelin–de Donder formula (Eqs. 4 and 6): with $A^f$ the affinity of the forward reaction and $A^r$ that of the reverse one,
#
# ```math
# v = \kappa\left(e^{A^f/RT} - e^{A^r/RT}\right),
# ```
#
# with a rate constant $\kappa$ (mol/s). For a reaction $\mathrm{A} \rightleftharpoons \mathrm{B}$ the affinities are the species' potentials, $A^f = \mu_a$ and $A^r = \mu_b$, and substituting the potentials gives first-order mass action, $v = \kappa(K_a x_a - K_b x_b) = k^+ x_a - k^- x_b$ with $k^+ = \kappa K_a$ and $k^- = \kappa K_b$ (Eqs. 7–9).
#
# Here the two reactions have rates $v_1 = \kappa_1(K_a x_a - K_b x_b)$ and $v_2 = \kappa_2(K_b x_b - K_c x_c)$. Reaction 1 consumes A and produces B, reaction 2 consumes B and produces C, so the mass balances are (Eq. 30)
#
# ```math
# \dot{x}_a = -v_1, \qquad \dot{x}_b = v_1 - v_2, \qquad \dot{x}_c = v_2.
# ```
#
# Their sum is zero, so the total amount $x_a + x_b + x_c$ is conserved. At equilibrium both flows vanish, $K_a x_a = K_b x_b = K_c x_c$, so the equilibrium constants of the two reactions are
#
# ```math
# K_1 = \frac{x_b}{x_a} = \frac{K_a}{K_b} = \frac{k_1^+}{k_1^-}, \qquad K_2 = \frac{x_c}{x_b} = \frac{K_b}{K_c} = \frac{k_2^+}{k_2^-},
# ```
#
# and each species ends with a share of the total proportional to $1/K_i$.
#
# See also P. J. Gawthrop and E. J. Crampin, [Energy-based analysis of biochemical cycles using bond graphs](https://doi.org/10.1098/rspa.2014.0459), *Proc. R. Soc. A* 470: 20140459, 2014 (open access; [arXiv:1406.2447](https://arxiv.org/abs/1406.2447), whose section and equation numbers this example cites), Fig. 1 for the bond graphs of this and other small networks.
#
# ### Bond Graph Model
#
# Each species is a C element, a store of matter whose effort is its chemical potential: it charges as a capacitor does when molar flow comes in. Each reaction is an Re element, a two-port that relates its flow to the affinities on its two bonds and dissipates power, $(A^f - A^r)\,v \ge 0$. A 0-junction (common potential) at each species lets several reactions share it. The bond graph is drawn from the model below with `bgplot`, once it is built.
#
# | Element | Stands for | Parameters |
# |---|---|---|
# | `Ce:A`, `Ce:B`, `Ce:C` | the pools of A, B, C | $K_a$, $K_b$, $K_c$, $RT$ |
# | `Re:r1`, `Re:r2` | the reactions A ⇌ B and B ⇌ C | $\kappa_1$, $\kappa_2$, $RT$ |
#
# ## Preliminary Steps
#
# In this example we first define the two chemical elements, which the library doesn't have, then build and simulate the bond graph model and check its equilibrium against the equilibrium constants. Finally, we simulate the mass-action equations derived above as a plain ModelingToolkit model and compare the results.
#
# We import the bond-graph toolkit to model the system and `OrdinaryDiffEqTsit5` to solve the resulting ODE. We also import ModelingToolkit's independent variable `t` and derivative `D`, for the custom elements and the reference model.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using Plots

# ## System Element Definitions
#
# ### Chemical species
#
# A species is a nonlinear compliance (Karnopp, Margolis & Rosenberg 2012, §3.1, pp. 40–42: a C relates its effort to its displacement $q = \int f\,dt$). Here the displacement is the amount $x$ and the effort the chemical potential $\mu = RT\ln(K x)$. `@oneport` adds the power port `power`, and the law is written in its bond variables, as for the built-in elements. We call the constructor `Ce` to keep it apart from the library's linear `Compliance`; Gawthrop & Crampin (2014) label it `C`.
@oneport function Ce(; name, K = 1.0, RT = 2577.0, x = 1.0)
    @parameters K = K           # (1/mol) Thermodynamic constant, K = exp(μ⁰/RT)/V
    @parameters RT = RT         # (J/mol) Gas constant times temperature
    @variables x(t) = x         # (mol) Amount of the species, with its initial value
    [D(x) ~ f, e ~ RT * log(K * x)]
end

# ### Reactions
#
# A reaction is a two-port: port 1 is the bond into it, carrying the forward affinity $A^f$, and port 2 the bond out of it, carrying the reverse affinity $A^r$. Both bonds carry the reaction flow $v$, which the Marcelin–de Donder formula sets. `@twoport` gives the bond variables `e1, f1, e2, f2`, each flow positive along its half-arrow:
@twoport function Re(; name, κ = 1.0, RT = 2577.0)
    @parameters κ = κ           # (mol/s) Reaction rate constant
    @parameters RT = RT         # (J/mol) Gas constant times temperature
    [f1 ~ κ * (exp(e1 / RT) - exp(e2 / RT)), f2 ~ f1]
end

# ### Elements
#
# Body temperature, $T = 310$ K, gives $RT$:
RT = 8.314 * 310.0              # (J/mol) R = 8.314 J/(mol⋅K)

# The species, each with its thermodynamic constant and initial amount. A chemical potential is $-\infty$ at a zero amount, so every species starts with some:
@named A = Ce(K = 1.0, RT = RT, x = 0.8)     # (1/mol, J/mol, mol)
@named B = Ce(K = 2.0, RT = RT, x = 0.1)
@named C = Ce(K = 0.5, RT = RT, x = 0.1)

# The two reactions:
@named r1 = Re(κ = 1.0, RT = RT)             # (mol/s) A ⇌ B
@named r2 = Re(κ = 0.5, RT = RT)             # (mol/s) B ⇌ C

# ## Model Construction
#
# Each species sits on a 0-junction, which gives it free ports for the reactions. Following Gawthrop & Crampin (2014, Fig. 1c and Eq. 10), every bond to a species points out of its junction into the species (a bare member `x`), so the species' flow is its rate of change, $\dot{x}_i$; this is also Karnopp's usual orientation for a C (§3.1, p. 42). Each reaction's bonds follow the forward reaction: from A's junction into `r1` (free port `:r1`, leaving), out of `r1` into B's junction (`:r1 => :in`, entering), and the same for `r2`. B's junction thus has one bond in and two out, and its flow balance is $v_1 - \dot{x}_b - v_2 = 0$, the mass balance of B.
@named ja = Junction0(A, :r1)
@named jb = Junction0(:r1 => :in, B, :r2)
@named jc = Junction0(:r2 => :in, C)

# A reaction's port 1 is its bond in and port 2 its bond out, so each junction's port towards a reaction connects to `p1` and the one from it to `p2`:
@named network = System([connect(ja.r1, r1.p1), connect(r1.p2, jb.r1),
                         connect(jb.r2, r2.p1), connect(r2.p2, jc.r2)], t;
                        systems = [ja, jb, jc, r1, r2])

# The library draws the bond graph from the model, as in Gawthrop & Crampin's Fig. 1c:
bgplot(network)

# ## Analysing the model
#
# `mtkcompile` reduces the model to its states: the amounts of the three species.
sys = mtkcompile(network)

# Print system states
unknowns(sys)

# Print system parameters
parameters(sys)

# Print the system's equations in terms of the states
full_equations(sys)

# ### Generated equations against the hand derivation
#
# `bg_equations` returns the same equations with the names we choose. We rename the amounts and parameters to the symbols of the introduction (inside a `let`, so the names `A`, `B` and `C` keep pointing at the elements); every element's $RT$ becomes the one $RT$:
let
    @variables x_a(t) x_b(t) x_c(t)
    @parameters K_a K_b K_c κ_1 κ_2 RT
    bg_equations(network; rename = Dict(
        ja.A.x => x_a, jb.B.x => x_b, jc.C.x => x_c,
        ja.A.K => K_a, jb.B.K => K_b, jc.C.K => K_c, r1.κ => κ_1, r2.κ => κ_2,
        ja.A.RT => RT, jb.B.RT => RT, jc.C.RT => RT, r1.RT => RT, r2.RT => RT))
end

# With one $RT$, the potential's $RT$ cancels the reaction's, and each exponential reads $e^{\ln(K_i x_i)} = K_i x_i$: the forward and reverse affinities of each reaction are the potentials of its reactant and product, set by the 0-junctions. The three equations are therefore $\dot{x}_a = -\kappa_1(K_a x_a - K_b x_b)$, $\dot{x}_b = \kappa_1(K_a x_a - K_b x_b) - \kappa_2(K_b x_b - K_c x_c)$ and $\dot{x}_c = \kappa_2(K_b x_b - K_c x_c)$, the mass-action equations of the introduction.
#
# ## System Simulations
#
# ### Bond graph model
#
# Define the simulation time, long enough to reach equilibrium
tspan = (0.0, 50.0)

# Generate an `ODEProblem`, with the initial amounts given to the elements:
prob = ODEProblem(sys, [], tspan)

# Use `solve` to simulate the system. Tight tolerances let us compare it with the reference model at the end.
sol = solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-10)
plot(sol, idxs = [ja.A.x, jb.B.x, jc.C.x], labels = ["x_a" "x_b" "x_c"],
     xlabel = "Time (s)", ylabel = "Amount (mol)")

# ### Equilibrium against the equilibrium constants
#
# The total amount is 1 mol, so the equilibrium amounts are $x_i = (1/K_i)/\sum_j (1/K_j)$, with ratios $K_1 = K_a/K_b$ and $K_2 = K_b/K_c$:
Ka, Kb, Kc = 1.0, 2.0, 0.5
x_eq = [1 / Ka, 1 / Kb, 1 / Kc] / (1 / Ka + 1 / Kb + 1 / Kc)

# The bond graph's final amounts, their ratios against $K_1$ and $K_2$, and the largest distance from the predicted equilibrium:
x_end = sol[[ja.A.x, jb.B.x, jc.C.x]][end]
(x_end[2] / x_end[1] - Ka / Kb, x_end[3] / x_end[2] - Kb / Kc, maximum(abs, x_end - x_eq))

# Every reaction dissipates, as the Re element should: the power into it, $(A^f - A^r)\,v$, is never negative along the trajectory (Gawthrop & Crampin 2014, §2.1):
minimum(sol[(r1.p1.e - r1.p2.e) * r1.p1.f]), minimum(sol[(r2.p1.e - r2.p2.e) * r2.p1.f])

# ### Mass-action model
#
# The reference model is the mass-action equations derived by hand, with the rate constants $k^\pm$ of each reaction, $k_1^+ = \kappa_1 K_a$, $k_1^- = \kappa_1 K_b$, $k_2^+ = \kappa_2 K_b$ and $k_2^- = \kappa_2 K_c$ (Eq. 9):
@variables xa(t) = 0.8 xb(t) = 0.1 xc(t) = 0.1      # (mol) Amounts of A, B, C
k1p, k1m = 1.0 * Ka, 1.0 * Kb                        # (1/s) Rate constants of A ⇌ B
k2p, k2m = 0.5 * Kb, 0.5 * Kc                        # (1/s) Rate constants of B ⇌ C

# The two reaction rates and the mass balances:
v1 = k1p * xa - k1m * xb
v2 = k2p * xb - k2m * xc
@named ref = System([D(xa) ~ -v1, D(xb) ~ v1 - v2, D(xc) ~ v2], t)

# We compile and solve it as the bond graph model, with the same time span and tolerances.
ref_sol = solve(ODEProblem(mtkcompile(ref), [], tspan), Tsit5(), reltol = 1e-10, abstol = 1e-10);

# ### Compare models
#
# The bond graph (solid) against the mass-action model (dashed):
plot(sol, idxs = [ja.A.x, jb.B.x, jc.C.x], labels = ["x_a, bond graph" "x_b, bond graph" "x_c, bond graph"],
     xlabel = "Time (s)", ylabel = "Amount (mol)")
plot!(ref_sol, idxs = [xa, xb, xc], labels = ["x_a, mass action" "x_b, mass action" "x_c, mass action"],
      linestyle = :dash)

# On a common time grid, the largest difference between the two models, for each species:
tgrid = range(tspan..., length = 1001)
[maximum(abs, sol(tgrid; idxs = bg).u - ref_sol(tgrid; idxs = ma).u)
 for (bg, ma) in [(ja.A.x, xa), (jb.B.x, xb), (jc.C.x, xc)]]

# The bond graph model reproduces mass-action kinetics to the solver's tolerance, and settles at the equilibrium the equilibrium constants predict.
