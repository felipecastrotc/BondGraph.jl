# # [Signals and Detectors](@id signals)
#
# Bond graphs carry information on *activated bonds*: signals that drive a source or a two-port's modulus, or that report a measured effort or flow, without carrying power (Karnopp, Margolis & Rosenberg, §3.2; Borutzky 2011, Ch. 3). Here the signals are MTK standard-library blocks:
#
# - `Se()`, `Sf()` without a value read their effort or flow from `input` (MSe, MSf);
# - `TF()`, `GY()` without a modulus read it from `input` (MTF, MGY);
# - `De()` and `Df()` put the effort or flow of their bond on `output`.

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
sol = solve(ODEProblem(mtkcompile(circuit), [], (0.0, 10.0)), Tsit5())
plot(sol, idxs = [rc.C.q], xlabel = "Time", ylabel = "Charge")

# ## Detectors
#
# A flow detector sits on a 1-junction and an effort detector on a 0-junction, where every bond shares the variable they read. Their bond carries no power, so they don't change the model. Here `Df` reads the circuit's current:
@named V = Se(2.0)
@named R = Resistance(R = 4.0)
@named ammeter = Df()
@named circuit = Junction1(V, [-1, R], [-1, ammeter])
sol = solve(ODEProblem(mtkcompile(circuit), [], (0.0, 1.0)), Tsit5())
sol[circuit.ammeter.output.u][end]   # E/R = 0.5

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
sol = solve(ODEProblem(mtkcompile(loop), [], (0.0, 10.0)), Rodas5P())
plot(sol, idxs = [out.C.q], xlabel = "Time", ylabel = "Displacement")
