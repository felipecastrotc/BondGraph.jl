# # [Pendulum Driven by a DC Motor](@id motor-pendulum)
#
# ## Introduction
#
# A DC motor swinging a pendulum is the smallest model of a robot joint: an electrical drive, a rotor, and a link whose geometry changes as it moves. Breedveld uses it to show how a bond graph is built from an iconic diagram and what causality tells the modeller (Breedveld 2011, in Borutzky (ed.), *Bond Graph Modelling of Engineering Systems*, Springer, §1.3.5, §1.4.4 and §1.5.9, pp. 19–35). The iconic diagram is his Fig. 1.6 (p. 19): a voltage source $u_s$ drives a current through the armature resistance $R_e$ and inductance $L$ into a motor with constant $K_m$; the rotor has inertia $J$ and bearing friction $b$; a rigid massless rod of length $l$ carries a point mass $m$ in gravity $g$.
#
# ### Motor and pendulum model
#
# As in the [DC motor](@ref dc-motor), the magnetic field is constant, so the motor is an ideal transducer with torque $K_m i$ and back EMF $K_m \omega$. The rod is rigid and massless, the mass a point, and the friction viscous. The angle $\varphi$ is measured from the horizontal axis $x$, and gravity acts along $-y$ (Fig. 1.6), so the mass sits at
#
# $x = l\cos\varphi, \qquad y = l\sin\varphi$
#
# (Eqs. 1.13–1.14, p. 20) and hangs at rest at $\varphi = -\pi/2$.
#
# Kirchhoff's voltage law around the armature loop gives
#
# $L\,\dfrac{di}{dt} + R_e\,i = u_s - K_m\,\omega.$
#
# The rotor and the mass turn together about the pivot, with total inertia $J + m l^2$. Gravity $-mg$ along $y$ acts at the lever arm $x = l\cos\varphi$, a torque $-m g l\cos\varphi$. Newton's 2nd law for the rotation gives
#
# $(J + m l^2)\,\dot{\omega} + b\,\omega + m g l\cos\varphi = K_m\,i, \qquad \dot{\varphi} = \omega.$
#
# Breedveld reaches the same model from the bond graph: the two inertias of the point mass add up to $m l^2$ on the rotor (Eqs. 1.22–1.23, p. 34), and gravity becomes a nonlinear compliance with torque $m g l\cos\varphi$ (Eq. 1.25, p. 34; Fig. 1.11b, p. 35). For small angles around $-\pi/2$ this compliance is linear, with parameter $1/(mgl)$, and the undamped pendulum oscillates at
#
# $f = \dfrac{\omega_0}{2\pi} = \dfrac{1}{2\pi}\sqrt{\dfrac{m g l}{J + m l^2}}$
#
# (Eq. 1.26, p. 35).
#
# ### Bond Graph Model
#
# Breedveld's bond graph (Fig. 1.7d, p. 20) keeps the point mass in the plane instead of reducing it by hand. Here is the model built below, drawn by `bgplot`:
#
# - `Se:us` is the voltage source $u_s$; `R:Re` and `I:L` the armature resistance $R_e$ and inductance $L$, on the 1-junction of the current $i$;
# - `GY:Km` is the motor, a gyrator with modulus $K_m$: $e_1 = K_m f_2$, $e_2 = K_m f_1$ (Breedveld Eqs. 1.6–1.7, p. 15; Karnopp, Margolis & Rosenberg 2012, Eq. 3.11);
# - `I:J` and `R:b` are the rotor inertia and the bearing friction, on the 1-junction of the angular velocity $\omega$;
# - `TF:mtfx` and `TF:mtfy` are the rigid rod: two transformers modulated by the angle, with moduli $-l\sin\varphi$ and $l\cos\varphi$ (Eq. 1.19, p. 21). With the transformer law $e_1 = n e_2$, $f_2 = n f_1$ (Eqs. 1.4–1.5, p. 15) they give the velocities of the mass, $\dot{x} = (-l\sin\varphi)\,\omega$ and $\dot{y} = (l\cos\varphi)\,\omega$ (Eqs. 1.15–1.16), and the torque of the rod's forces, $T_\varphi = (-l\sin\varphi)\,F_x + (l\cos\varphi)\,F_y$ (Eq. 1.18);
# - `I:mx` and `I:my` are the point mass in the two directions, and `Se:mg` the weight $-mg$ on the 1-junction of $\dot{y}$.
#
# The angle that modulates the transformers is the integral of $\omega$, a signal: Breedveld draws it as an integrator fed from the $\omega$ junction (Fig. 1.7d). Here a flow detector `Df:ω` reads $\omega$ and an integrator block outputs $\varphi$. Signals are not drawn by `bgplot`.
#
# ## Preliminary Steps
#
# In this example we build the bond graph of Fig. 1.7d and let ModelingToolkit reduce it, then compare it with the reduced pendulum derived above and check the small-angle frequency of Eq. 1.26.
#
# We import the bond-graph toolkit, the integrator block from the MTK standard library, and the solvers: the compiled model keeps algebraic equations (see below), so it needs the DAE-capable `Rodas5P`, and `OrdinaryDiffEqNonlinearSolve` for its initial values. `Tsit5` solves the reference ODE.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Blocks as Blocks
using OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using Plots

# ## System Element Definitions
#
# The electrical side: the voltage source, the armature resistance and inductance.
@named us = Se(1.0)                 # (V) Source voltage uₛ
@named Re = Resistance(R = 1.0)     # (Ω) Armature resistance Rₑ
@named L = Inertia(I = 0.01)        # (H) Armature inductance L

# The motor constant $K_m$, in N⋅m/A = V⋅s/rad, is the gyrator's modulus:
@named Km = GY(r = 0.5)

# The rotor: its inertia and the bearing friction.
@named J = Inertia(I = 0.01)        # (kg⋅m²) Rotor inertia J
@named b = Resistance(R = 0.05)     # (N⋅m⋅s/rad) Bearing friction b

# The pendulum's mass, rod length and gravity:
m = 0.5                             # (kg) Point mass m
l = 0.3                             # (m) Rod length l
g = 9.81                            # (m/s²) Gravitational acceleration g

# The point mass is one inertia per direction. Its momenta are fixed by the rotor speed through the rigid rod, so they are not independent states (below) and get no initial value, `p = nothing`; an initial value would over-determine the model's initialisation.
@named mx = Inertia(I = m, p = nothing)   # (kg) Point mass, x direction
@named my = Inertia(I = m, p = nothing)   # (kg) Point mass, y direction

# Gravity is a constant force $-mg$ along $y$ (Fig. 1.7d labels it `Se: −mg`):
@named mg = Se(-m * g)              # (N) Weight

# The angle: a flow detector reads $\omega$ on its 1-junction, and an integrator turns it into $\varphi$. Its initial value is the starting angle, the pendulum hanging at rest:
@named ω = Df()
@named φ = Blocks.Integrator(x = -π / 2)

# The rod's two transformers read the integrator's output in their moduli, so they are built with `@named` (see [Transformers and Gyrators](@ref two-ports)):
@named mtfx = TF(m = -l * sin(φ.output.u))
@named mtfy = TF(m = l * cos(φ.output.u))

# ## Model Construction
#
# The bonds follow Fig. 1.7d. The source's bond points into the current junction and the bonds of $R_e$ and $L$ out of it, into the elements, as in the [DC motor](@ref dc-motor); the bond to the motor leaves the junction (`[-1, :gy]`):
@named ji = Junction1(us, [-1, Re], [-1, L], [-1, :gy])

# The motor's bond enters the $\omega$ junction (`:gy`); the bonds of $J$, $b$ and the detector point out of it, and so do the two bonds to the rod's transformers (`[-1, :x]`, `[-1, :y]`): power flows from the motor to the mass.
@named jω = Junction1(:gy, [-1, J], [-1, b], [-1, ω], [-1, :x], [-1, :y])

# The $y$ transformer's bond enters the junction of $\dot{y}$ (`:tf`), and so does the weight's; the bond of `my` points out of it. In the $x$ direction there is nothing else on the velocity $\dot{x}$, so Breedveld bonds the transformer straight into `mx` (Fig. 1.7d), as here.
@named jy = Junction1(:tf, mg, [-1, my])

# Each transformer's port 1 is the bond into it and port 2 the bond out of it. The last `connect` is the signal from the detector to the integrator.
@named model = System([connect(ji.gy, Km.p1), connect(Km.p2, jω.gy),
                       connect(jω.x, mtfx.p1), connect(mtfx.p2, mx.power),
                       connect(jω.y, mtfy.p1), connect(mtfy.p2, jy.tf),
                       connect(jω.ω.output, φ.input)], t;
                      systems = [ji, Km, jω, mtfx, mx, mtfy, jy, φ])

# The bond graph, drawn by the library:
bgplot(model)

# ## Analysing the model
#
# Breedveld's causality assignment leaves both inertias of the point mass in derivative causality (Fig. 1.8, p. 26): the rigid rod fixes their velocities by the rotor's, so the three inertias are not independent (§1.5.9, p. 33). The library assigns no causality; `mtkcompile` finds the same dependency and reduces the index of the model.
sys = mtkcompile(model)

# Print system states
unknowns(sys)

# Print system parameters
parameters(sys)

# Print the system's equations in terms of the states
full_equations(sys)

# ### Generated equations against the hand derivation
#
# `bg_equations` gives the same equations with names of our choice. ModelingToolkit kept the angle $\varphi$ (the integrator's output), the $x$-momentum of the mass $p_x$, the rotor speed $\omega$, the flux linkage $\lambda = L\,i$, the rod's force on the $\dot{y}$ junction $F_y$ (the effort of the $y$ transformer's port 2) and $\dot{\omega}$, a *dummy derivative* added by the index reduction. We rename them, and the parameters, to the symbols of the introduction; the integrator's gain is 1 (inside a `let`, so the names `L`, `J`, `b`, `m` keep pointing at the elements and values):
let
    @variables φ(t) p_x(t) ω(t) λ(t) F_y(t) ω̇(t)
    @parameters R_e L K_m J b m
    ω_t = only(filter(v -> endswith(string(v), "fˍt(t)"), unknowns(sys)))   # the dummy derivative
    bg_equations(model; rename = Dict(
        model.φ.output.u => φ, model.mx.p => p_x, model.jω.J.power.f => ω,
        model.ji.L.p => λ, model.jy.core.p1.e => F_y, ω_t => ω̇,
        model.φ.k => 1, model.ji.Re.R => R_e, model.ji.L.I => L, model.Km.r => K_m,
        model.jω.J.I => J, model.jω.b.R => b, model.mx.I => m, model.jy.my.I => m))
end

# The rod length and the weight appear as their values, $l = 0.3$ and $mg = 4.905$, and the source as $u_s = 1$. Read against the hand derivation:
#
# - $\dot{\varphi} = \omega$, the integrator, and $\dot{\lambda} = u_s - R_e\,i - K_m\,\omega$ with $i = \lambda/L$, Kirchhoff's law of the introduction.
# - $0 = -p_x/m - l\sin\varphi\,\omega$ is the $x$ transformer's flow law, $\dot{x} = (-l\sin\varphi)\,\omega$ (Eq. 1.15), and the equation for $\dot{p}_x$ its derivative, $\dot{p}_x = m\,\ddot{x}$. These are how MTK recovers $\omega$ from the state $p_x$.
# - The $\dot{y}$ junction gives $\dot{p}_y = F_y - mg$, with $p_y = m\,\dot{y} = m\,l\cos\varphi\,\omega$ differentiated: $(F_y - mg)/m = l\cos\varphi\,\dot{\omega} - l\sin\varphi\,\omega^2$.
# - The $\omega$ junction is the torque balance $K_m i = J\dot{\omega} + b\,\omega + (l\cos\varphi)\,F_y + (-l\sin\varphi)\,\dot{p}_x$, the rod's torque of Eq. 1.18 with the mass's forces $F_x = \dot{p}_x$ and $F_y$.
#
# Substituting $\dot{p}_x$ and $F_y$ into the torque balance, the centripetal terms in $\omega^2$ cancel (Breedveld's gyristors, Eq. 1.24) and the mass terms add up to $m l^2 \dot{\omega}$ (Eq. 1.22), leaving $(J + m l^2)\,\dot{\omega} + b\,\omega + m g l\cos\varphi = K_m i$: the hand-derived equation.
#
# **Note:** ModelingToolkit chooses the states once, at compile time. With $p_x$ as a state, $\omega$ is recovered from $p_x = -m l\sin\varphi\,\omega$, which fails where $\sin\varphi = 0$, the rod horizontal. The runs below stay in the lower half, $-\pi < \varphi < 0$.
#
# ## System Simulations
#
# ### Bond graph model
#
# Define the simulation time
tspan = (0.0, 10.0)

# Generate an `ODEProblem`. All states start from the defaults: no current, the pendulum at rest at $\varphi = -\pi/2$.
prob = ODEProblem(sys, [], tspan)

# Use `solve` to simulate the system. The source's voltage drives a torque $K_m u_s / R_e = 0.5$ N⋅m in steady state, which swings the pendulum up to where gravity balances it, $m g l\cos\varphi = 0.5$.
sol = solve(prob, Rodas5P(), reltol = 1e-10, abstol = 1e-10)
plot(sol, idxs = [φ.output.u, jω.J.power.f, ji.L.power.f],
     labels = ["φ (rad)" "ω (rad/s)" "i (A)"], xlabel = "Time (s)")

# ### Reduced pendulum model
#
# The hand-derived equations of the introduction, as a plain MTK `System`, with the same parameter values and initial state:
@variables i(t) = 0.0           # Armature current
@variables ϕ(t) = -π / 2        # Pendulum angle
@variables ω̃(t) = 0.0           # Angular velocity

Reᵥ = 1.0                       # (Ω) Armature resistance
Lᵥ = 0.01                       # (H) Armature inductance
Uᵥ = 1.0                        # (V) Source voltage
Kᵥ = 0.5                        # (N⋅m/A) Motor constant
Jᵥ = 0.01                       # (kg⋅m²) Rotor inertia
bᵥ = 0.05                       # (N⋅m⋅s/rad) Bearing friction

eqs_t = [
    Lᵥ * D(i) ~ Uᵥ - Reᵥ * i - Kᵥ * ω̃,                           # Armature loop
    (Jᵥ + m * l^2) * D(ω̃) ~ Kᵥ * i - bᵥ * ω̃ - m * g * l * cos(ϕ),  # Rotation
    D(ϕ) ~ ω̃,                                                    # Angle
]

# Then, we create the `System` and compile it.
@named sys_t = System(eqs_t, t)
sys_t = mtkcompile(sys_t)

# The reduced model is an ODE, solved with `Tsit5` and the same tolerances.
sol_t = solve(ODEProblem(sys_t, [], tspan), Tsit5(), reltol = 1e-10, abstol = 1e-10);

# ### Compare models
#
# Pendulum angle
plot(sol, idxs = [φ.output.u], label = "Bond graph", ylabel = "Angle (rad)")
plot!(sol_t, idxs = [ϕ], label = "Reduced model", xlabel = "Time (s)", linestyle = :dash)

# Armature current
plot(sol, idxs = [ji.L.power.f], label = "Bond graph", ylabel = "Current (A)")
plot!(sol_t, idxs = [i], label = "Reduced model", xlabel = "Time (s)", linestyle = :dash)

# On a common time grid, the largest difference between the two models, in angle, angular velocity and current:
tgrid = range(tspan..., length = 1001)
(maximum(abs, sol(tgrid; idxs = φ.output.u).u - sol_t(tgrid; idxs = ϕ).u),
 maximum(abs, sol(tgrid; idxs = jω.J.power.f).u - sol_t(tgrid; idxs = ω̃).u),
 maximum(abs, sol(tgrid; idxs = ji.L.power.f).u - sol_t(tgrid; idxs = i).u))

# ### Small-angle frequency
#
# Without the motor's coupling and the friction ($K_m = 0$, $b = 0$), the pendulum released at a small angle $\psi_0$ from $-\pi/2$ swings as $\psi_0\cos\omega_0 t$, with $\omega_0$ from Eq. 1.26. The same compiled model takes the new parameter values and starting angle:
ψ₀ = 1e-3
prob_0 = ODEProblem(sys, [Km.r => 0.0, jω.b.R => 0.0, φ.x => -π / 2 + ψ₀], tspan)
sol_0 = solve(prob_0, Rodas5P(), reltol = 1e-12, abstol = 1e-12)

# The frequency of Eq. 1.26, in Hz:
ω₀ = sqrt(m * g * l / (Jᵥ + m * l^2))
ω₀ / 2π

# The simulated deviation $\psi = \varphi + \pi/2$, scaled by $\psi_0$, against $\cos\omega_0 t$:
ψ = (sol_0(tgrid; idxs = φ.output.u).u .+ π / 2) ./ ψ₀
plot(tgrid, ψ, label = "Bond graph", xlabel = "Time (s)", ylabel = "ψ / ψ₀")
plot!(τ -> cos(ω₀ * τ), 0, 10, label = "cos ω₀t", linestyle = :dash)

# Their largest difference. It comes from the pendulum's nonlinearity, not the solver: a pendulum's frequency drops with amplitude, $\omega \approx \omega_0 (1 - \psi_0^2/16)$ (the standard large-amplitude correction, not from Breedveld), a phase lag of $\psi_0^2 \omega_0 t / 16 \approx 3 \cdot 10^{-6}$ after 10 s.
maximum(abs, ψ - cos.(ω₀ .* tgrid))

# The bond graph of Fig. 1.7d, with the point mass in derivative causality, reproduces the reduced pendulum to the solver's tolerance, and its small oscillations have the frequency of Eq. 1.26.
