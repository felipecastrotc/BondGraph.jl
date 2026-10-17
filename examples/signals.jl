# # [Signals and Detectors](@id signals)
#
# Bond graphs carry information on *activated bonds*: signals that drive a source or a two-port's modulus, or that report a measured effort or flow, without carrying power (Karnopp, Margolis & Rosenberg, §3.2; Borutzky 2011, Ch. 3). Here the signals are MTK standard-library blocks:
#
# - `Se()`, `Sf()` without a value read their effort or flow from `input` (MSe, MSf);
# - `TF()`, `GY()` without a modulus read it from `input` (MTF, MGY);
# - `De()` and `Df()` put the effort or flow of their bond on `output`.
#
# Each model below is checked against its closed-form solution.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t
import ModelingToolkitStandardLibrary.Blocks as Blocks
using OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using Plots

# ## A source driven by a block
#
# A sine voltage from `Blocks.Sine` drives a series RC circuit: the source's bond points into the 1-junction, the R and C bonds point out of it.
@named V = Se()
@named R = Resistance(R = 1.0)
@named C = Compliance(C = 1.0)
@named rc = Junction1(V, [-1, R], [-1, C])
@named sine = Blocks.Sine(frequency = 0.5)
@named circuit = System([connect(sine.output, rc.V.input)], t; systems = [rc, sine])
sol = solve(ODEProblem(mtkcompile(circuit), [], (0.0, 10.0)), Tsit5(), reltol = 1e-10, abstol = 1e-10)

# By hand: the 1-junction shares the current $\dot{q}$, and its efforts sum to zero with `+` for the source's bond and `-` for the R and C bonds, $V(t) - R\,\dot{q} - q/C = 0$. With $R = C = 1$, $V = \sin\omega t$ and $\omega = 2\pi \cdot 0.5 = \pi$, the charge from $q(0) = 0$ is
#
# $q(t) = \dfrac{\sin\omega t - \omega\cos\omega t + \omega\,e^{-t}}{1 + \omega^2}.$
ω = π
q_rc(τ) = (sin(ω * τ) - ω * cos(ω * τ) + ω * exp(-τ)) / (1 + ω^2)
plot(sol, idxs = [rc.C.q], label = "Bond graph", xlabel = "Time", ylabel = "Charge")
plot!(q_rc, 0, 10, label = "Closed form", linestyle = :dash)

# The largest difference on a common time grid:
tgrid = 0:0.01:10
maximum(abs, sol(tgrid; idxs = rc.C.q).u - q_rc.(tgrid))

# ## Detectors
#
# A flow detector sits on a 1-junction and an effort detector on a 0-junction, where every bond shares the variable they read. Their bond carries no power, so they don't change the model. Here `Df` reads the circuit's current, which is $E/R = 2/4 = 0.5$ at all times:
@named V = Se(2.0)
@named R = Resistance(R = 4.0)
@named ammeter = Df()
@named circuit = Junction1(V, [-1, R], [-1, ammeter])
sol = solve(ODEProblem(mtkcompile(circuit), [], (0.0, 1.0)), Tsit5())
sol[circuit.ammeter.output.u][end]

# ## Closing a loop
#
# A detector, blocks and a signal modulus close a feedback loop entirely in MTK. A source pushes flow through a gyrator into a compliance; an effort detector reads the compliance's effort, and the blocks set the gyrator's modulus to `r = 2 + 0.5 e`:
@named C = Compliance(C = 1.0)
@named de = De()
@named out = Junction0(:in, [-1, C], [-1, de])
@named se = Se(2.0)
@named src = Junction1(se, [-1, :out])
@named gy = GY()
@named gain = Blocks.Gain(k = 0.5)
@named offset = Blocks.Constant(k = 2.0)
@named add = Blocks.Add()
@named loop = System([connect(src.out, gy.p1), connect(gy.p2, out.in),
                      connect(out.de.output, gain.input), connect(gain.output, add.input1),
                      connect(offset.output, add.input2), connect(add.output, gy.input)], t;
                     systems = [src, gy, out, gain, offset, add])

# A modulus read from a signal is an unknown of the model, which MTK won't divide by, so the compiled model keeps an algebraic equation. Solve it with a DAE-capable solver such as `Rodas5P`; `OrdinaryDiffEqNonlinearSolve` provides the nonlinear solve for its initial values.
sol = solve(ODEProblem(mtkcompile(loop), [], (0.0, 10.0)), Rodas5P(), reltol = 1e-10, abstol = 1e-10)

# By hand: the source's 1-junction passes the effort $e_1 = 2$ to the gyrator's port 1. The gyrator law $e_1 = r f_2$ (Karnopp, Margolis & Rosenberg 2012, Eq. 3.11) sets the flow out of port 2, $f_2 = 2/r$, and the 0-junction passes it to the compliance (the detector takes no flow), so $\dot{q} = 2/r$. The detector reads $e = q/C = q$, so $r = 2 + q/2$ and
#
# $\dot{q} = \dfrac{4}{4 + q}, \qquad q(0) = 0 \quad\Rightarrow\quad q(t) = \sqrt{16 + 8t} - 4.$
q_loop(τ) = sqrt(16 + 8τ) - 4
plot(sol, idxs = [out.C.q], label = "Bond graph", xlabel = "Time", ylabel = "Displacement")
plot!(q_loop, 0, 10, label = "Closed form", linestyle = :dash)

# The largest difference on a common time grid:
maximum(abs, sol(tgrid; idxs = out.C.q).u - q_loop.(tgrid))
