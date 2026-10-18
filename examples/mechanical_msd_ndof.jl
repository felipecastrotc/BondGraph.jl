# # [MSD system with N-DoF](@id msd-ndof)
#
# ## Introduction
#
# One of the advantages of a scripting approach to modeling is that a model can be built algorithmically. To showcase this with the bond-graph toolkit, we consider a mass–spring–damper (MSD) chain with $N$ degrees of freedom, built by repeating one bond pattern in a comprehension.
#
# ### N-DoF MSD model
#
# The chain has $N$ masses $m_1, \dots, m_N$ in a line, each moving along the line with displacement $x_i$ from its rest position (positive in the same direction for all masses). The first mass is attached to the ground by a spring $k_1$ and a damper $c_1$ in parallel. Each pair of neighbouring masses $i-1$ and $i$ is coupled by a spring $k_i$ and a damper $c_i$ in parallel, $i = 2, \dots, N$. The springs and dampers are linear and massless, and no external force acts on the chain.
#
# Treat the ground as a mass $0$ that does not move, $x_0 = 0$. Spring–damper pair $i$ joins mass $i-1$ to mass $i$; when $x_{i-1} - x_i > 0$ it is compressed, so it pushes mass $i$ forward and mass $i-1$ back with the force
#
# $F_i = k_i\,(x_{i-1} - x_i) + c_i\,(\dot{x}_{i-1} - \dot{x}_i), \qquad i = 1, \dots, N.$
#
# Mass $i$ is pushed by pair $i$ behind it and pulled back by pair $i+1$ ahead of it (there is no pair $N+1$, so $F_{N+1} = 0$). Newton's 2nd law then gives, for every mass,
#
# $m_i\,\ddot{x}_i = F_i - F_{i+1}, \qquad i = 1, \dots, N.$
#
# For $N = 4$, written out:
#
# $m_1\,\ddot{x}_1 = -k_1\,x_1 - c_1\,\dot{x}_1 - k_2\,(x_1 - x_2) - c_2\,(\dot{x}_1 - \dot{x}_2),$
#
# $m_2\,\ddot{x}_2 = k_2\,(x_1 - x_2) + c_2\,(\dot{x}_1 - \dot{x}_2) - k_3\,(x_2 - x_3) - c_3\,(\dot{x}_2 - \dot{x}_3),$
#
# $m_3\,\ddot{x}_3 = k_3\,(x_2 - x_3) + c_3\,(\dot{x}_2 - \dot{x}_3) - k_4\,(x_3 - x_4) - c_4\,(\dot{x}_3 - \dot{x}_4),$
#
# $m_4\,\ddot{x}_4 = k_4\,(x_3 - x_4) + c_4\,(\dot{x}_3 - \dot{x}_4).$
#
# All masses, stiffnesses and damping coefficients are one. The chain starts at rest with the spring between masses 1 and 2 stretched by one unit: $x_1 = 0$ and $x_2 = x_3 = x_4 = 1$, so the other springs are at their rest length.
#
# ### Bond Graph Model
#
# Karnopp, Margolis & Rosenberg (2012) model a mass as an inertia `I:m` with momentum $p = m\,\dot{x}$ and $\dot{p} = F$, a spring as a compliance `C:1/k` with displacement $q$ and effort $e = q/C = k\,q$, and a dashpot as a resistance `R:c` (§3.1, pp. 38–47). A mechanical 1-junction is a force balance at a common velocity, and a mechanical 0-junction expresses the geometric compatibility of velocities at a common force (§3.3, pp. 58, 61).
#
# The bond graph has one 1-junction per velocity:
#
# - the velocity $\dot{x}_1$ of mass 1, with `I:m₁`, `C:1/k₁` and `R:c₁` on it;
# - the velocity $\dot{x}_i$ of mass $i = 2, \dots, N$, with `I:mᵢ` on it;
# - the relative velocity of spring–damper pair $i = 2, \dots, N$, with `C:1/kᵢ` and `R:cᵢ` on it.
#
# Between the junctions of masses $i-1$ and $i$ sits a 0-junction, which shares the force of pair $i$ among the three 1-junctions around it and relates their velocities. In short, the repeating pattern is `1 (mass i−1) → 0 → 1 (mass i)`, with a third bond `0 → 1 (pair i)`. The [2-DoF](@ref msd-2dof) example draws the same pattern for two masses, with a force on the second.
#
# ## Preliminary Steps
#
# In this example we build and simulate the bond graph model of the chain for $N = 4$. Then we simulate the hand-derived equations above as a plain ModelingToolkit model and compare the two.
#
# We import the bond-graph toolkit to model the system and `OrdinaryDiffEqTsit5` to solve the resulting ODE. We also import ModelingToolkit's independent variable `t` and derivative `D` for the reference model at the end.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using Plots

# Set the number of degrees of freedom. The construction below needs $N \geq 2$.
n = 4

# ## System Element Definitions
#
# The elements of mass 1: its inertia, the spring to ground (a compliance $C = 1/k_1$) and the damper to ground:
@named m = Inertia(I = 1.0)       # (kg) Mass
@named s = Compliance(C = 1.0)    # (m/N) Spring compliance, 1/k
@named d = Resistance(R = 1.0)    # (N⋅s/m) Damping coefficient

# The masses $2, \dots, N$ reuse `m`, and the coupling dampers reuse `d`. The coupling springs are created per pair below, because the second one starts with an initial displacement.

# ## Model Construction
#
# ### Mass 1
#
# Mass 1 is a 1-junction: its inertia, spring and damper share its velocity $\dot{x}_1$. Following the usual orientation (Karnopp, Margolis & Rosenberg 2012, §3.1, pp. 39, 42), the bonds of the passive elements point out of the junction into the elements (a bare member `x`). A free port named `out`, also pointing out of the junction (`:out`), will carry the force of the first coupling. The junction then states $0 = -\dot{p}_1 - k_1 q_1 - c_1\,\dot{x}_1 - e_{\text{out}}$, where the compliance's displacement $q_1 = x_1$ is the stretch of the spring to ground.
@named m1 = Junction1(m, d, s, :out)

# ### Masses 2 to N
#
# Each further mass is a 1-junction with a free port `in`, pointing into the junction from the coupling behind it (`:in => :in`), and its inertia `m`. Every mass but the last also has a port `out`, pointing out towards the next coupling. So mass $i$ states $\dot{p}_i = e_{\text{in}} - e_{\text{out}}$: the force from behind minus the force ahead.
masses = [Junction1(:in => :in, m, (i < n ? [:out] : [])...; name = Symbol(:m, i))
          for i in 2:n]

# ### Spring–damper pairs
#
# Pair $i$ is a 1-junction whose free port `j0` points into it from the 0-junction (`:j0 => :in`), with its compliance and resistance bonds pointing out into the elements. It states $e_{\text{j0}} = q_i/C_i + c_i\,f$: the force carried to the 0-junction is the spring force plus the damper force, at the pair's relative velocity $f$.
#
# The compliance of pair 2 starts at $q = -1$. The sign follows from the bond directions, as the 0-junction below shows: they make $q_i = x_{i-1} - x_i$, which is negative when the spring is stretched. The initial stretch of one unit is therefore $q_2 = -1$.
sds = [Junction1(Compliance(C = 1.0, q = i == 2 ? -1.0 : 0.0, name = :s), d, :j0 => :in;
                 name = Symbol(:sd, i))
       for i in 2:n]

# ### Couplings
#
# Each coupling is a 0-junction joining three free ports: `out` of mass $i-1$, `in` of mass $i$ and `j0` of pair $i$. A plain `connect` of free ports forms that 0-junction: the efforts are shared and the signed flows sum to zero, each flow counted positive when its bond points into the 0-junction (Karnopp, Margolis & Rosenberg 2012, Eqs. 3.16–3.17, pp. 58–59). The signs come from the half-arrows chosen above, not from the order of the `connect` arguments.
#
# The bond from mass $i-1$ points into the 0-junction and the bonds to mass $i$ and to pair $i$ point out of it. Flow continuity gives $\dot{x}_{i-1} = \dot{x}_i + f$, so the pair's velocity is $f = \dot{x}_{i-1} - \dot{x}_i$ and its displacement is $q_i = x_{i-1} - x_i$: positive when the spring is compressed. The common effort is then the force $F_i$ of the hand derivation, which enters mass $i$ with a plus sign (`in`) and mass $i-1$ with a minus sign (`out`), as in $m_i\,\ddot{x}_i = F_i - F_{i+1}$.
prev = [m1; masses[1:end-1]]
cons = [connect(prev[k].out, masses[k].in, sds[k].j0) for k in 1:n-1]
@named model = System(cons, t; systems = [m1; masses; sds])

# ## Analysing the model
#
# `mtkcompile` reduces the model to its states: the momenta of the $N$ inertias and the displacements of the $N$ compliances.
sys = mtkcompile(model)

# Print system states
unknowns(sys)

# Print system parameters
parameters(sys)

# Print the system's equations in terms of the states
full_equations(sys)

# ### Generated equations against the hand derivation
#
# `bg_equations` returns the same equations with the names we choose. Here, for $N = 4$, we rename the momenta to $p_i$, the displacements to $q_i$ and the parameters to $m_i$, $C_i = 1/k_i$ and $c_i$ (inside a `let`, so these helper names stay local):
let
    @variables p₁(t) p₂(t) p₃(t) p₄(t) q₁(t) q₂(t) q₃(t) q₄(t)
    @parameters m₁ m₂ m₃ m₄ C₁ C₂ C₃ C₄ c₁ c₂ c₃ c₄
    J = [m1; masses]
    p, q, M, C, c = [p₁, p₂, p₃, p₄], [q₁, q₂, q₃, q₄], [m₁, m₂, m₃, m₄], [C₁, C₂, C₃, C₄], [c₁, c₂, c₃, c₄]
    S = [m1; sds]
    rename = Dict([[J[i].m.p => p[i] for i in 1:n]; [J[i].m.I => M[i] for i in 1:n];
                   [S[i].s.q => q[i] for i in 1:n]; [S[i].s.C => C[i] for i in 1:n];
                   [S[i].d.R => c[i] for i in 1:n]])
    bg_equations(model; rename)
end

# MTK writes the equations without collecting terms, so some appear twice with opposite signs, and $\dot{x}_i = p_i/m_i$ shows up as a quotient. With $p_i = m_i\,\dot{x}_i$, $q_1 = x_1$, $q_i = x_{i-1} - x_i$ and $q_i/C_i = k_i\,q_i$:
#
# - the equations for $\dot{q}_1$ and $\dot{q}_i$ restate $\dot{q}_1 = \dot{x}_1$ and $\dot{q}_i = \dot{x}_{i-1} - \dot{x}_i$;
# - the equation for $\dot{p}_1$ is $m_1\,\ddot{x}_1 = -k_1 x_1 - c_1\,\dot{x}_1 - F_2$;
# - the equation for $\dot{p}_i$, $1 < i < N$, is $m_i\,\ddot{x}_i = F_i - F_{i+1}$, and the one for $\dot{p}_N$ is $m_N\,\ddot{x}_N = F_N$,
#
# where, for $i \geq 2$, $F_i = q_i/C_i + c_i\,(p_{i-1}/m_{i-1} - p_i/m_i) = k_i\,(x_{i-1} - x_i) + c_i\,(\dot{x}_{i-1} - \dot{x}_i)$: the governing equations of the introduction.

# ## System Simulations
#
# ### Bond graph model
#
# Define the simulation time
tspan = (0.0, 10.0)

# Generate an `ODEProblem` as for any ModelingToolkit model. The initial values are those set on the elements: zero momenta, and the initial displacement of the second compliance.
prob = ODEProblem(sys, [], tspan)

# Use `solve` to simulate the system. Tight tolerances let us compare it with the reference model at the end.
sol = solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-10)

# The velocity of each mass is the flow of its inertia's bond:
velocities = [J.m.f for J in [m1; masses]]
plot(sol, idxs = velocities, xlabel = "Time (s)", ylabel = "Velocity (m/s)",
     label = permutedims(["ẋ$i" for i in 1:n]))

# ### N-DoF MSD model
#
# The reference model is the hand-derived $m_i\,\ddot{x}_i = F_i - F_{i+1}$ as a plain ModelingToolkit `System`, with the displacements $x_i$ and velocities $v_i = \dot{x}_i$ as states. First, the parameters (all one, as in the bond graph model):
mᵥ = ones(n)                      # (kg) Masses
kᵥ = ones(n)                      # (N/m) Spring stiffnesses
cᵥ = ones(n)                      # (N⋅s/m) Damping coefficients

# Then the states, with the initial stretch of the spring between masses 1 and 2:
@variables x(t)[1:n] = [0.0; ones(n - 1)]   # Displacements
@variables v(t)[1:n] = zeros(n)             # Velocities

# The force $F_i$ of spring–damper pair $i$, with the ground as $x_0 = \dot{x}_0 = 0$ and no pair beyond the last mass:
F(i) = i > n ? 0 :
       i == 1 ? -kᵥ[1] * x[1] - cᵥ[1] * v[1] :
       kᵥ[i] * (x[i-1] - x[i]) + cᵥ[i] * (v[i-1] - v[i])

# Newton's 2nd law for every mass:
eqs_t = [[D(x[i]) ~ v[i] for i in 1:n];
         [mᵥ[i] * D(v[i]) ~ F(i) - F(i + 1) for i in 1:n]]

# Then, we use the equations to create a `System` and compile it with `mtkcompile`, as for the bond graph model.
@named sys_t = System(eqs_t, t)
sys_t = mtkcompile(sys_t)

# Define the `ODEProblem` with the same time span and solve it with the same tolerances.
prob_t = ODEProblem(sys_t, [], tspan)
sol_t = solve(prob_t, Tsit5(), reltol = 1e-10, abstol = 1e-10)
plot(sol_t, idxs = collect(v), xlabel = "Time (s)", ylabel = "Velocity (m/s)",
     label = permutedims(["v$i" for i in 1:n]))

# ### Compare models
#
# The velocities of the bond graph model (solid) and of the reference model (dashed):
plot(sol, idxs = velocities, label = permutedims(["Bond graph ẋ$i" for i in 1:n]))
plot!(sol_t, idxs = collect(v), linestyle = :dash, label = permutedims(["Model v$i" for i in 1:n]),
      xlabel = "Time (s)", ylabel = "Velocity (m/s)")

# The displacements of the compliances against those of the springs in the reference model, $q_1 = x_1$ and $q_i = x_{i-1} - x_i$:
stretch = [m1.s.q; [S.s.q for S in sds]]
stretch_t = [x[1]; [x[i-1] - x[i] for i in 2:n]]
plot(sol, idxs = stretch, label = permutedims(["Bond graph q$i" for i in 1:n]))
plot!(sol_t, idxs = stretch_t, linestyle = :dash, label = permutedims(["Model q$i" for i in 1:n]),
      xlabel = "Time (s)", ylabel = "Displacement (m)")

# On a common time grid, the largest difference between the two models over all velocities and spring displacements:
tgrid = range(tspan..., length = 1001)
maximum(maximum(abs, sol(tgrid; idxs = a).u - sol_t(tgrid; idxs = b).u)
        for (a, b) in zip([velocities; stretch], [collect(v); stretch_t]))

# The bond graph model, built by repeating one pattern, reproduces the hand-derived dynamics of the chain to the solver's tolerance.
