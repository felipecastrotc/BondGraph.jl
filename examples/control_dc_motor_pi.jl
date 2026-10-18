# # [DC Motor Speed Control](@id control-dc-motor-pi)
#
# ## Introduction
#
# A DC motor left to itself settles at whatever speed its voltage, its resistance and its friction allow. To make it turn at a commanded speed, a controller measures the speed, compares it with the reference and adjusts the armature voltage. This example closes that loop around the motor of the [DC Motor](@ref dc-motor) example: a flow detector measures the shaft speed, an MTK standard-library PI controller turns the speed error into a voltage, and a modulated effort source applies it. The bond graph carries the power; the loop is made of signals, which carry none.
#
# The problem is the one of the University of Michigan *Control Tutorials for MATLAB and Simulink* (CTMS), [DC Motor Speed: System Modeling](https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling) and [DC Motor Speed: PID Controller Design](https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=ControlPID), whose physical sketch is
#
# ![DC motor](https://ctms.engin.umich.edu/CTMS/Content/MotorSpeed/System/Modeling/figures/motor.png)
#
# ### Plant
#
# CTMS assumes a rigid rotor and shaft, viscous friction, and a constant magnetic field, so that the motor torque is $T = K i$ and the back EMF is $e = K\dot{\theta}$, with one constant $K$ for both in SI units. Newton's second law and Kirchhoff's voltage law give
#
# $J\ddot{\theta} + b\,\dot{\theta} = K\,i,$
#
# $L\frac{di}{dt} + R\,i = V - K\,\dot{\theta},$
#
# and, eliminating $I(s)$, the transfer function from the armature voltage to the shaft speed $\omega = \dot{\theta}$,
#
# $P(s) = \dfrac{\Omega(s)}{V(s)} = \dfrac{K}{(Js + b)(Ls + R) + K^2}.$
#
# The parameters are CTMS's: $J = 0.01$ kg⋅m², $b = 0.1$ N⋅m⋅s, $K = 0.01$ (V⋅s/rad = N⋅m/A), $R = 1$ Ω, $L = 0.5$ H. They are the same values as in the [DC Motor](@ref dc-motor) example; what changes is the input: there a constant 12 V and an applied torque, here a voltage set by the controller and no load torque.
#
# ### Requirements and controller
#
# Without control, 1 V gives 0.1 rad/s in steady state. For a 1 rad/s step in the speed reference, CTMS asks for a settling time under 2 s, an overshoot under 5 % and a steady-state error under 1 %. Its controller is $C(s) = K_p + K_i/s + K_d s$ acting on the error $\omega_r - \omega$, with unity feedback; its PID page ends with $K_p = 100$, $K_i = 200$, $K_d = 10$.
#
# This example keeps the integral term, which removes the steady-state error, and leaves out the derivative: a PI controller, $V = K_p(\omega_r - \omega) + K_i\int(\omega_r - \omega)\,dt$. Its gains are not CTMS's; we choose them here. The plant's denominator is $JL\,s^2 + (JR + bL)\,s + bR + K^2 = 0.005\,s^2 + 0.06\,s + 0.1001$, with poles near $-2$ and $-10$. Putting the PI zero $-K_i/K_p$ on the slow pole cancels it, and the loop becomes $\dfrac{2K_p}{s(s + 10)}$, whose closed loop $s^2 + 10\,s + 2K_p$ is a second-order system with $\omega_n^2 = 2K_p$ and $\zeta = 5/\omega_n$. With $K_p = 20$ and $K_i = 40$, $\zeta \approx 0.79$: about 2 % overshoot and a settling time near $4/(\zeta\omega_n) = 0.8$ s.
#
# The closed loop has three states, the speed, the current and the integral of the error $z$:
#
# $J\dot{\omega} = K\,i - b\,\omega,$
#
# $L\frac{di}{dt} = V - R\,i - K\,\omega,$
#
# $\dot{z} = \omega_r - \omega, \qquad V = K_p(\omega_r - \omega) + K_i\,z.$
#
# ### Bond Graph Model
#
# The plant is the DC motor's bond graph without the torque source: an electrical 1-junction with the source $V$, the resistance $R$ and the inertia $L$; a mechanical 1-junction with the inertia $J$ and the resistance $b$; and the gyrator, with modulus $K$, between them. The source is now modulated (MSe): its effort comes from a signal. A flow detector Df on the mechanical junction reads the shared flow, the shaft speed, and outputs it as a signal. Its bond carries no power ($e = 0$), so it doesn't change the plant. The feedback path, from Df through the error and the PI controller to the source, consists of activated bonds (Karnopp, Margolis & Rosenberg 2012, §3.2).
#
# The drawing is made from the model with `bgplot`, once the model is built.
#
# ## Preliminary Steps
#
# We build the closed loop from bond-graph elements and standard-library blocks, simulate a 1 rad/s step in the reference, and compare it with the closed-loop equations above written as a plain MTK model.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Blocks as Blocks
using OrdinaryDiffEqTsit5
using Plots

# ## System Element Definitions
#
# CTMS's motor parameters, and the PI gains chosen above. Both models below take their numbers from here.
p = (J = 0.01,     # (kg⋅m²) rotor moment of inertia
     b = 0.1,      # (N⋅m⋅s) viscous friction
     K = 0.01,     # (V⋅s/rad = N⋅m/A) motor constant
     R = 1.0,      # (Ω) armature resistance
     L = 0.5,      # (H) armature inductance
     Kp = 20.0,    # (V⋅s/rad) proportional gain
     Ki = 40.0,    # (V/rad) integral gain
     ωr = 1.0)     # (rad/s) speed reference

# The electrical side: the armature inductance and resistance, and the armature voltage as a modulated effort source, `Se()` without a value, which reads its effort from its `input`:
@named L = Inertia(I = p.L)
@named R = Resistance(R = p.R)
@named V = Se()

# The mechanical side: the rotor inertia, its friction, and the speed sensor, a flow detector whose `output` is the flow of its bond:
@named J = Inertia(I = p.J)
@named b = Resistance(R = p.b)
@named tach = Df()

# The motor's gyrator, with the motor constant as its modulus:
@named gy = GY(r = p.K)

# The control blocks: the reference, the error $\omega_r - \omega$ (`Feedback` outputs `input1 - input2`), and the PI controller. `Blocks.PI` computes $k\,(e + \frac{1}{T}\int e\,dt)$, so $k = K_p$ and $T = K_p/K_i$.
@named ref = Blocks.Constant(k = p.ωr)
@named err = Blocks.Feedback()
@named pi = Blocks.PI(k = p.Kp, T = p.Kp / p.Ki)

# ## Model Construction
#
# The junctions follow the DC motor example. On each 1-junction the source's bond points in, away from the source, and the bonds of R, I and the detector point out into the element (Karnopp, Margolis & Rosenberg 2012, §3.1, pp. 39 and 49). The electrical bond leaves its junction towards the gyrator and the mechanical bond enters its junction from it, following the power from the source to the shaft (Fig. 2.3d, pp. 25–26). So the back EMF enters the electrical effort sum with a minus sign and the motor torque enters the mechanical one with a plus sign.
@named je = Junction1(V => :pout, R, L, :gy)
@named jm = Junction1(:gy => :in, J, b, tach)

# The bonds go through the gyrator, port 1 in and port 2 out; the signals go from the detector to the error's second input, from the reference to its first, then through the controller to the source.
@named loop = System([connect(je.gy, gy.p1), connect(gy.p2, jm.gy),
                      connect(ref.output, err.input1), connect(jm.tach.output, err.input2),
                      connect(err.output, pi.err_input), connect(pi.ctr_output, je.V.input)], t;
                     systems = [je, gy, jm, ref, err, pi])

# The bond graph, drawn from the model. `bgplot` draws the bonds only: the signal path from `tach` through the blocks to `V` is not drawn, and the modulated source is labelled `Se`.
bgplot(loop)

# ## Analysing the model
#
# `mtkcompile` keeps three states: the two inertias' momenta and the controller's integrator.
sys = mtkcompile(loop)

# Print system states
unknowns(sys)

# Print system parameters
parameters(sys)

# ### Generated equations against the hand derivation
#
# `bgequations` names the states and parameters as we choose. The inertias' states are their momenta, the flux linkage $\lambda = L\,i$ and the angular momentum $h = J\omega$; the controller's state is its integrator's output, $x = z/T$. We rename them, and the parameters, to the symbols of the introduction, inside a `let` so that `L`, `R`, `J` and `b` keep naming the elements:
let
    @variables λ(t) h(t) x(t)
    @parameters L R J b K Kp T ωᵣ
    bgequations(loop; rename = Dict(je.L.p => λ, jm.J.p => h, pi.addPI.input2.u => x,
                                     je.L.I => L, je.R.R => R, jm.J.I => J, jm.b.R => b,
                                     gy.r => K, pi.k => Kp, pi.T => T, ref.k => ωᵣ,
                                     pi.addPI.k1 => 1, pi.addPI.k2 => 1))
end

# With $\lambda = L\,i$ and $h = J\omega$, these read $L\frac{di}{dt} = V - R\,i - K\omega$ and $J\dot{\omega} = K\,i - b\,\omega$, the plant's equations, where the source's effort is the controller's output $V = K_p\,(\omega_r - \omega + x)$ (the weights of `Blocks.PI`'s internal `Add`, both 1, are substituted), and the integrator gives $T\dot{x} = \omega_r - \omega$. With $x = z/T$ and $K_p/T = K_i$, this is the closed loop of the introduction.

# ## System Simulations
#
# ### Bond graph model
#
# Everything starts at rest, with the reference at 1 rad/s: a step at $t = 0$.
tspan = (0.0, 3.0)
prob = ODEProblem(sys, [], tspan)
sol = solve(prob, Tsit5(), reltol = 1e-10, abstol = 1e-10)
plot(sol, idxs = [jm.J.f, je.L.f], labels = ["ω (rad/s)" "i (A)"], xlabel = "Time (s)")

# ### Closed-loop reference model
#
# The closed-loop equations of the introduction, as a plain MTK model with CTMS's parameters and the same gains:
@variables ω(t) = 0.0 i(t) = 0.0 z(t) = 0.0 Vₐ(t)
eqs_cl = [
    p.J * D(ω) ~ p.K * i - p.b * ω,
    p.L * D(i) ~ Vₐ - p.R * i - p.K * ω,
    D(z) ~ p.ωr - ω,
    Vₐ ~ p.Kp * (p.ωr - ω) + p.Ki * z,
]
@named sys_cl = System(eqs_cl, t)
sys_cl = mtkcompile(sys_cl)
sol_cl = solve(ODEProblem(sys_cl, [], tspan), Tsit5(), reltol = 1e-10, abstol = 1e-10);

# ### Compare models
#
# Shaft speed:
plot(sol.t, sol[jm.J.f], label = "Bond graph", ylabel = "Speed (rad/s)")
plot!(sol_cl.t, sol_cl[ω], label = "Closed-loop equations", xlabel = "Time (s)", linestyle = :dash)

# Armature voltage, the controller's output:
plot(sol.t, sol[pi.ctr_output.u], label = "Bond graph", ylabel = "Voltage (V)")
plot!(sol_cl.t, sol_cl[Vₐ], label = "Closed-loop equations", xlabel = "Time (s)", linestyle = :dash)

# On a common time grid, the largest difference between the two models in speed, current and voltage:
tgrid = range(tspan..., length = 1001)
(maximum(abs, sol(tgrid; idxs = jm.J.f).u - sol_cl(tgrid; idxs = ω).u),
 maximum(abs, sol(tgrid; idxs = je.L.f).u - sol_cl(tgrid; idxs = i).u),
 maximum(abs, sol(tgrid; idxs = pi.ctr_output.u).u - sol_cl(tgrid; idxs = Vₐ).u))

# ### Against CTMS's requirements
#
# The overshoot (%), the 2 % settling time (s) and the steady-state error (%) at the end of the run, for a 1 rad/s step:
ωg = sol(tgrid; idxs = jm.J.f).u
overshoot = 100 * (maximum(ωg) - p.ωr) / p.ωr
settling = tgrid[findlast(w -> abs(w - p.ωr) > 0.02 * p.ωr, ωg) + 1]
steady_error = 100 * abs(p.ωr - ωg[end]) / p.ωr
(overshoot, settling, steady_error)

# All three meet CTMS's requirements (under 5 %, under 2 s, under 1 %). The bond graph with a detector, blocks and a modulated source reproduces the hand-derived closed loop to the solver's tolerance.
