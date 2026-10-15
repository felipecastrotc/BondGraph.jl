# Research 11: How do power ports couple with ModelingToolkitStandardLibrary components?

Ticket: `.scratch/mtk-upgrade/issues/11-coupling-with-mtk-stdlib.md` · Date: 2026-10-15

## Summary

1. **A direct `connect(power, <MTKSL connector>)` is impossible in every domain, and MTK reports it as an error, not a silent mis-pairing.** MTK requires the connected connectors' unknowns to have the same *names*. `Power(e, f)` against `Pin(v, i)`, `MechanicalPort(v, f)`, `Flange(phi, tau)`, `Flange(s, f)`, `HydraulicPort(p, dm)` or `RealInput(u)` throws "Different types of connectors" **[ran]**. Silent positional mis-pairing only happens when the names are the same set in a different order (research 03, F6).
2. **Renaming the variables would not help either.** In the mechanical domains MTKSL's potential is position or velocity and its flow variable is force or torque. That is the *mobility* analogy: the bond-graph effort is MTKSL's *through* variable. In the electrical domain, `Pin.v` is an absolute node potential, not a port voltage **[reasoned; ran for a renamed (v, i) port: it works, but only as a grounded one-port]**.
3. **A tiny adapter element works in every domain.** It has one `Power` port and the MTKSL connector(s), takes 2–3 equations, needs no MTK internals, and conserves power. The MTKSL network attached to it *is* the bond-graph element at the end of the bond: MTKSL's through-variable is negated; the across-variable maps to the effort (electrical, hydraulic) or, through its derivative, to the flow (rotational, translational).
4. **Verified on MTK 11.45.4 / MTKSL 2.29.8 against reference models** **[ran]**:
   - RC: BG vs pure MTKSL, ≤ 3.5e-10.
   - DC motor: BG vs MTKSL load, the reverse direction, and BG vs a pure MTKSL `EMF` motor, ≤ 5.7e-10.
   - Translational: velocity- and position-based, ≤ 7.7e-11.
   - Hydraulic: ≤ 4e-11 relative.
   - `RealInput`-modulated GY: bitwise identical, or 7.7e-13 in closed loop through blocks.
5. **Blocks.** A source or a GY/TF ratio driven by a block is a component with a `RealInput`, and a sensor is a zero-effort or zero-flow element with a `RealOutput` **[ran]**. The existing "ratio is an expression" route also works (`@named gy = GY(g = blk.output.u)`), but only because `@named` wraps symbolic keyword arguments in `ParentScope`. Without `@named` it fails with an unbalanced-system error **[ran]**.
6. **Recommendation:**
   - Add five named per-domain adapters: `ElectricalAdapter`, `TranslationalAdapter`, `TranslationalPositionAdapter`, `RotationalAdapter`, `HydraulicAdapter`. Name their ports as MTKSL does (`p`/`n`, `flange`, `port`).
   - Add a `RealInput` form for `Se`/`Sf`/`GY`/`TF` and `EffortSensor`/`FlowSensor`.
   - Make MTKSL a regular dependency. It is MIT-licensed, depends only on ModelingToolkitBase, and adds ~0.2 s to load.
   - Watch the export clashes `Power` (`Blocks.Power`) and `Mass`/`Spring`/`Damper` (Translational/Rotational).

Environment: Julia 1.10.12, ModelingToolkit 11.45.4, ModelingToolkitBase 1.77.4, ModelingToolkitStandardLibrary 2.29.8 (latest release; tag commit `319811b`), OrdinaryDiffEqRosenbrock (`Rodas5P`, abstol = reltol = 1e-10 unless noted). Bond-graph elements come from the accepted prototype `prototype/connection_mechanism/BGProto.jl`, which was `include`d and not modified.

Source abbreviations:
- **MSL** = `https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/`
- **MTKB** = `https://github.com/SciML/ModelingToolkit.jl/blob/v11.45.4/lib/ModelingToolkitBase/src/`

The cited lines were checked against both tags.

---

## 1. Why a direct `connect` cannot work

| # | Fact | Evidence |
|---|---|---|
| D1 | Connection expansion takes the first connector's unknowns and errors if any other connector has a different count **or an unknown whose name is not in that set** (`any(!in(sys1_dvs_set), unknown_vars)`). It then checks that the connect type (Equality/Flow) is equal per position. Pairing is by position *after* the name check. | [MTKB systems/connectors.jl#L649-L652](https://github.com/SciML/ModelingToolkit.jl/blob/v11.45.4/lib/ModelingToolkitBase/src/systems/connectors.jl#L649-L652), [#L670-L673](https://github.com/SciML/ModelingToolkit.jl/blob/v11.45.4/lib/ModelingToolkitBase/src/systems/connectors.jl#L670-L673), error text [#L370-L372](https://github.com/SciML/ModelingToolkit.jl/blob/v11.45.4/lib/ModelingToolkitBase/src/systems/connectors.jl#L370-L372) |
| D2 | **[ran]** `connect(bg.power, X)` for X = `Electrical.Pin`, `Translational.MechanicalPort`, `Rotational.Flange`, `TranslationalPosition.Flange`, `HydraulicPort`, `Blocks.RealInput` → every one errors with `Different types of connectors are in one connection statement`. No MTKSL connector has unknowns named `{e, f}` (all `@connector`s in `src/` grepped: also `MagneticPort(V_m, Phi)` and `HeatPort(T, Q_flow)`), so a silent mis-pair with `Power` cannot happen today. | `coupling/probe_direct.jl` |
| D3 | **[ran]** A port with the *same* names as `Pin` (`v`, `i[Flow]`) does connect. It then behaves as a **grounded** one-port: its effort is the node potential against ground (series 2 Ω + "BG" 5 Ω at 3 V gives i = 3/7 exactly). This is what a renamed `Power` would give in the electrical domain. A floating element needs two pins. | `probe_direct.jl` (VIPower part) |
| D4 | In mechanics, MTKSL is position- or velocity-across and force-through: `Flange(phi, tau[Flow])`, `Flange(s, f[Flow])`, `MechanicalPort(v, f[Flow])`. In bond-graph terms the effort (force, torque) is MTKSL's *flow* variable, so no renaming or reordering can make the two connectors the same "type". A renamed `Power` would equate force with velocity. **[reasoned]** | [MSL Rotational/utils.jl#L10-L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Rotational/utils.jl#L10-L16); [TranslationalPosition/utils.jl#L10-L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/TranslationalPosition/utils.jl#L10-L16); [Translational/utils.jl#L10-L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Translational/utils.jl#L10-L16) |
| D5 | MTKSL's sign convention is Modelica's: every through-variable is positive **into** the component, and a through-input source writes the opposite sign (`flange.f ~ -f`). This is the same convention as `Power.f`, which is the main reason the adapters stay trivial. | [MSL docs/src/connectors/sign_convention.md#L1-L30](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/docs/src/connectors/sign_convention.md#L1-L30); `Pin` comment "Current flowing into the pin" [Electrical/utils.jl#L13](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Electrical/utils.jl#L13) |

## 2. The adapter convention

An adapter is an ordinary bond-graph **element** (ADR 0001: element bonds point into the element). The MTKSL subnetwork on its MTKSL side *is* that element:

- `power.f` is the flow into the subnetwork, and `power.e` is the effort across it.
- MTKSL's through-variable at the adapter's connector is "into the adapter". It is the negative of what flows into the MTKSL subnetwork, so it is **negated**.
- The across-variable maps to the effort in the impedance-analogy domains (electrical, hydraulic). In the mobility-analogy domains (mechanics) it maps, through its rate, to the flow.

Check: power from BG into the adapter, `e·f`, equals power out through the MTKSL connector, `−(across·through)`, in every case **[reasoned]**. All runs below confirm it.

**Transparency check [reasoned, then ran]:**
- An MTKSL `Resistor` across `p`/`n` gives exactly `e = R f`.
- An MTKSL `Voltage` source across `p`/`n` behaves exactly like `Se`: `e = V`, and the flow into it is negative while it delivers power.
- An MTKSL `Inertia` on `flange` gives `J dω/dt = e`, `f = ω`.

So **in a `Junction1`, give the adapter the sign you would give the element it stands in for**: +1 where an `Se` would sit (the MTKSL side supplies power, as in examples 1a/5a), and −1 where an R/C/I would sit (the MTKSL side is a load, as in examples 2a/3). On a 0-junction (`connect`) no sign is needed. Wiring the adapter straight to a GY port also works (2b).

| Domain | MTKSL connector (MSL v2.29.8) | Effort / flow mapping | Adapter equations |
|---|---|---|---|
| Electrical | `Pin(v, i[Flow])`, i into pin. [Electrical/utils.jl#L10-L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Electrical/utils.jl#L10-L16); two-pin `OnePort` [#L33-L54](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Electrical/utils.jl#L33-L54) | e = voltage p→n, f = current through | `power.e ~ p.v - n.v`, `power.f ~ -p.i`, `0 ~ p.i + n.i` |
| Translational (velocity) | `MechanicalPort(v, f[Flow])`. [Translational/utils.jl#L10-L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Translational/utils.jl#L10-L16); `Mass`: `flange.f ~ f`, `D(v) ~ f/m` [components.jl#L76-L100](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Translational/components.jl#L76-L100) | e = force, f = velocity | `power.e ~ -flange.f`, `power.f ~ flange.v` |
| Translational (position) | `TranslationalPosition.Flange(s, f[Flow])` [utils.jl#L10-L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/TranslationalPosition/utils.jl#L10-L16); `TranslationalModelica.Flange(s, f[Flow])` has the same names and types [utils.jl#L10-L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/TranslationalModelica/utils.jl#L10-L16) | e = force, f = ds/dt | `power.e ~ -flange.f`, `D(flange.s) ~ power.f` |
| Rotational | `Flange(phi, tau[Flow])` [Rotational/utils.jl#L10-L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Rotational/utils.jl#L10-L16); `Inertia`: `J*a ~ flange_a.tau + flange_b.tau` [components.jl#L53-L80](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Rotational/components.jl#L53-L80) | e = torque, f = dφ/dt | `power.e ~ -flange.tau`, `D(flange.phi) ~ power.f` |
| Hydraulic | `HydraulicPort(p, dm[Flow])` plus 7 fluid **parameters** (ρ, β, μ, …) set by connecting a `HydraulicFluid` [IsothermalCompressible/utils.jl#L16-L33](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Hydraulic/IsothermalCompressible/utils.jl#L16-L33), [#L50-L74](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Hydraulic/IsothermalCompressible/utils.jl#L50-L74) | e = gauge pressure, f = volume flow at the reference density | `power.e ~ port.p`, `power.f ~ -port.dm / port.ρ` |

Notes per domain:

- **Electrical.** Use two pins, not one, because `Pin.v` is absolute (D3). A grounded use is `connect(adapter.n, ground.g)`. An MTKSL electrical node (a `connect` of pins) is a 0-junction **[reasoned]**.
- **Translational, velocity-based.** No derivative, so the adapter adds no state. An MTKSL translational node (common `v`, Σf = 0) is a mechanical 1-junction **[reasoned]**.
- **Position-based (rotational, `TranslationalPosition`, `TranslationalModelica`).** The adapter integrates the flow into an absolute angle or position. That adds one state the bond graph doesn't need, and **an initial value is required**, exactly as for a pure-MTKSL `Inertia`. Without `phi`, `ODEProblem` fails with "Initialization incomplete" **[ran: `probe_init.jl`]**. Don't give the adapter a hard default. When the MTKSL side also sets the value (`Mass(s = 0)`), MTK warns that initialization is over-determined **[ran: `probe_tmodelica.jl`]**. Take a `phi = nothing` / `s = nothing` keyword, like MTKSL components do. A moving support (`use_support`) is skipped; add it when a model needs the reaction on a housing.
- **Hydraulic.**
  - MTKSL's through-variable is **mass** flow. `p·dm` is not power, so the BG flow is taken as the volume flow `dm/ρ` at the reference density `port.ρ`. With MTKSL's compressible (Tait) density this is exact for the adapter's own equations but is an approximation of "volume flow" at high pressure **[reasoned]** ([`liquid_density`, utils.jl](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Hydraulic/IsothermalCompressible/utils.jl#L136-L140)).
  - **Gotcha [ran]:** `connect(fluid, src.port, adapter.port)` in a single statement leaves the system unbalanced. Connect the fluid on its own, as the MTKSL tests do: `connect(fluid, src.port), connect(src.port, adapter.port)` ([test/isothermal_compressible.jl#L37](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/test/isothermal_compressible.jl#L37)).
  - Hydraulic states are volumes (q ~ 1e-6 m³), so `abstol` must be ~1e-18. At 1e-10 the comparison was off by 2e-5 relative, purely from the tolerance **[ran]**.

## 3. Blocks (causal signals)

`RealInput(u)` / `RealOutput(u)` are scalar connectors with `input = true` / `output = true` metadata ([Blocks/utils.jl#L12-L33](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Blocks/utils.jl#L12-L33), [#L72-L93](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Blocks/utils.jl#L72-L93)). MTKSL's own block-driven sources expose such an input: `Voltage` has `V`, `Pressure` has `p`, `MassFlow` has `dm` ([Electrical/Analog/sources.jl#L16-L36](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Electrical/Analog/sources.jl#L16-L36), [Hydraulic/.../sources.jl#L11-L27](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Hydraulic/IsothermalCompressible/sources.jl#L11-L27)). The bond-graph equivalents are one-liners:

- **Driven source:** `SeIn` with `power.e ~ u.u`. An `Sf` version would be `power.f ~ -u.u`.
- **Modulated two-port:** `mGY` with a `RealInput` `g`: `p1.e ~ -g.u*p2.f`, `p2.e ~ g.u*p1.f`. `mTF` is analogous.
- **Sensors:** `FlowSensor` is a zero-effort element (`e ~ 0`, `y.u ~ f`) that sits on a 1-junction. `EffortSensor` is the dual (`f ~ 0`, `y.u ~ e`) on a 0-junction **[reasoned]**.

Expression route: the prototype already allows `GY(g = <any expression>)`. Referencing a block's output (`@named gy = GY(g = sn.output.u)`, with `sn` a sibling subsystem) works, because `@named` rewrites every symbolic keyword argument to `ParentScope` ([MTKB systems/abstractsystem.jl#L2856-L2868](https://github.com/SciML/ModelingToolkit.jl/blob/v11.45.4/lib/ModelingToolkitBase/src/systems/abstractsystem.jl#L2856-L2868), [#L2986-L2997](https://github.com/SciML/ModelingToolkit.jl/blob/v11.45.4/lib/ModelingToolkitBase/src/systems/abstractsystem.jl#L2986-L2997)). **Calling `GY(g = sn.output.u, name = :gy)` without `@named` namespaces the variable as `gy₊sn₊output₊u` and `mtkcompile` throws `ExtraVariablesSystemException`** **[ran: `probe_scope2.jl`, `ex4_mgy.jl` (4b'')]**. It fails loudly, not silently, but it is a trap. ADR 0001's "modulation is just an expression" should note it, or `@oneport`/two-port constructors should apply `ParentScope` themselves (`ModelingToolkit.default_to_parentscope` is internal: **UNVERIFIED** as API). Because `@named` does this wrapping, the prototype's `GY(g = 0.01 + 0.001 * jm.J.power.f)` works.

## 4. Results (all [ran], scripts in Appendix B)

| # | Model | Reference | max \|Δ\| |
|---|---|---|---|
| 1a | MTKSL `Voltage` ← `Blocks.Sine`, `ElectricalAdapter`, BG `Junction1(:src, [-1,R], [-1,C])` | pure MTKSL series RC (`Resistor`, `Capacitor`) | v_C 3.2e-10, i 1.6e-10 |
| 1b | MTKSL `Voltage` + `Resistor`, BG C only (adapter node→ground) | same | v_C 3.5e-10, i 1.7e-10 |
| 1c | BG RC driven by `SeIn` ← `Blocks.Sine` (no adapter) | same | v_C 3.2e-10 |
| 2a | BG DC motor (V, R, L, GY, load τ) → `RotationalAdapter` → MTKSL `Inertia` + `Damper` + `Fixed` | pure BG motor (prototype) | ω 5.7e-10, i 1.6e-10 |
| 2b | MTKSL `Voltage`+`Resistor`+`Inductor` → `ElectricalAdapter` → BG GY + BG mechanics | pure BG motor | ω 2.9e-10, i 1.5e-10 |
| 2c | pure MTKSL motor (`EMF` k = 0.01, `ConstantTorque`) | pure BG motor | ω 5.7e-10, i 1.4e-10 (independent sign check: BG GY ≡ MTKSL EMF) |
| 3a | BG m–k–d (forced by `SeIn` ← `Sine`) → `TranslationalAdapter` → MTKSL `Spring` → `Mass` | pure MTKSL Translational | v1 3.1e-11, v2 5.3e-11, x1 1.1e-11 |
| 3b | same through `TranslationalPositionAdapter` → `TranslationalPosition.Spring`/`Mass` | same | v1 3.8e-11, v2 7.7e-11 |
| 3c | `Se(1)` → `TranslationalPositionAdapter` → `TranslationalModelica.Mass(m=2)` | analytic v = t/2 | 3.9e-14 |
| 4a | `Blocks.Constant(0.01)` → `mGY.g` | `GY(g = 0.01)` | 0.0 (identical) |
| 4b | `Blocks.Sine` → `mGY.g` | `GY(g = 0.01 + 0.005 sin(πt))` (≠ unmodulated by 0.57) | 0.0 |
| 4b' | `@named gy = GY(g = sn.output.u)` | same | 0.0; without `@named`: error |
| 4c | ω → `FlowSensor` → `Gain(0.001)` → `Add(+0.01)` → `mGY.g` (closed loop) | prototype `GY(g = 0.01 + 0.001 ω)` | 7.7e-13 |
| 5a | MTKSL `Pressure` ← `Constant` + `HydraulicFluid` → `HydraulicAdapter` → BG series R + C | BG `Se` + R + C | 0.0; `dm = −ρQ` to 2.8e-17 |
| 5b | MTKSL `MassFlow` ← `Sine` → `HydraulicAdapter` → BG R ∥ C | BG `Sf(dm/ρ)` ∥ R ∥ C | 2.0e-8 Pa on 489 Pa (4e-11 rel., abstol 1e-18) |

State counts stay minimal. 1a/1b keep one state (`C.q`), 5a one, and 4c two. Position-based adapters add the absolute angle or position (2a: `inertia.phi`; 3b: `tpa.flange.s`), which is harmless to the physics. The only warnings came from redundant user `u0` in the MTKSL reference (3, 3c), not from the adapters.

## 5. Recommended API

**Named per-domain adapters, not a generic one.** A generic `Adapter(connector; effort = …, flow = …, differentiate = …)` would push the analogy choice and the signs (the only subtle part) onto the user, and it would not be shorter: each adapter is three lines. The per-domain names mirror MTKSL's module names, the port names mirror MTKSL's own (`p`/`n`, `flange`, `port`), and every adapter has a `power` port. Using one is then just `connect`, which is as close to MTK as possible:

```julia
@named ea = ElectricalAdapter()
connect(src.p, ea.p), connect(src.n, ea.n, gnd.g), connect(ea.power, rc.src)
```

```julia
# Bond-graph side: `power::Power`, an element (flow into it). MTKSL side: MTKSL connectors.
function ElectricalAdapter(; name)                       # Electrical.Pin p, n
    @named power = Power(); @named p = Electrical.Pin(); @named n = Electrical.Pin()
    System([power.e ~ p.v - n.v, power.f ~ -p.i, 0 ~ p.i + n.i], t; name, systems = [power, p, n])
end
function TranslationalAdapter(; name)                    # Translational.MechanicalPort (v, f)
    @named power = Power(); @named flange = Translational.MechanicalPort()
    System([power.e ~ -flange.f, power.f ~ flange.v], t; name, systems = [power, flange])
end
function TranslationalPositionAdapter(; name)            # TranslationalPosition/-Modelica Flange (s, f)
    @named power = Power(); @named flange = TranslationalPosition.Flange()
    System([power.e ~ -flange.f, D(flange.s) ~ power.f], t; name, systems = [power, flange])
end
function RotationalAdapter(; name)                       # Rotational.Flange (phi, tau), fixed support
    @named power = Power(); @named flange = Rotational.Flange()
    System([power.e ~ -flange.tau, D(flange.phi) ~ power.f], t; name, systems = [power, flange])
end
function HydraulicAdapter(; name)                        # HydraulicPort (p, dm); f = volume flow
    @named power = Power(); @named port = IsothermalCompressible.HydraulicPort()
    System([power.e ~ port.p, power.f ~ -port.dm / port.ρ], t; name, systems = [power, port])
end
```

The tested code is exactly this (Appendix A). Still to add before release: a `phi`/`s` initial-value keyword on the two position adapters (passed to the flange variable, default `nothing`); `@component`; and docstrings that state the "adapter = the element it replaces" sign rule.

**Blocks:**
- Sources and two-ports take either a value/expression, as today, or no value, in which case they expose `input::RealInput`. `Se()` gives `power.e ~ input.u`; `GY()` gives the ratio `input.u`.
- Add `EffortSensor` and `FlowSensor` with `output::RealOutput`.
- The `input`/`output` names follow MTKSL's SISO blocks. MTKSL's sources name the port after the quantity (`V`, `tau`, `dm`), but the bond-graph `e`/`f` would collide visually with `power.e`. This is a naming choice for the maintainer.

**Dependency.** Make MTKSL a regular `[deps]` entry rather than a package extension:
- **[ran]** Load cost after `using ModelingToolkit`: 0.19–0.20 s.
- **[read]** Its deps are only ModelingToolkitBase, Symbolics, SymbolicUtils, DiffEqBase, ChainRulesCore, PreallocationTools, PrecompileTools and LinearAlgebra, and its license is MIT ([Project.toml](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/Project.toml), [LICENSE](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/LICENSE)). It therefore does not worsen the deferred AGPL question from research 02.
- **[reasoned]** An extension cannot add new exported names to the parent module, so it would need stub functions in BondGraph, plus `RealInput` duplicated or deferred.

**Export clashes [read; effect reasoned from Julia's `using` rules, not run].**
- MTKSL `Blocks` exports `Power` ([Blocks.jl#L16](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Blocks/Blocks.jl#L16)).
- `Translational`/`TranslationalPosition` export `Mass, Spring, Damper`, and `Rotational` exports `Spring, Damper` ([Rotational.jl#L15](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/blob/v2.29.8/src/Mechanical/Rotational/Rotational.jl#L15)).

With `using BondGraph, ModelingToolkitStandardLibrary.Blocks`, an unqualified `Power` is ambiguous and errors on use. Options:
- Rename the connector `PowerPort`, which is the glossary term "Power port" anyway.
- Rename the mechanical one-ports to their glossary role names (`I`/`C`/`R`, or `Inertia`/`Compliance`/`Resistance`). Note that `Inertia` itself clashes with `Rotational.Inertia`.
- Document qualified imports (`import ModelingToolkitStandardLibrary as MSL`). The test scripts do this.

## 6. Open points for the maintainer

- Adopt the "adapter = the element it replaces" sign rule. The alternative, an adapter that acts as a source to the bond graph, is the same equations with (e, f) → (−e, −f). It conserves power too, but it would invert the junction-sign rule.
- Resolve the `Power`/`Mass`/`Spring`/`Damper` export clash (§5).
- Decide whether constructors should apply `ParentScope` to expression keyword arguments themselves, so that `GY(g = expr, name = :gy)` works without `@named` (§3).
- Not covered: Thermal `HeatPort(T, Q_flow)`, where `T·Q_flow` is not power (a pseudo-bond graph, so an adapter is ill-defined without entropy flow); Magnetic `MagneticPort(V_m, Phi)` (effort–flow-like, reasoned to work like electrical, not run); `MultiBody2D`.

---

## Appendix A: `coupling/adapters.jl` (shared by all scripts)

Run each script with:

```sh
julia --project=<env with ModelingToolkit 11.45.4, ModelingToolkitStandardLibrary 2.29.8,
                 OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock> coupling/<script>.jl
```

The env used was `/private/tmp/claude-501/-Users-fctc-me-doutorado-dev-bg-bg-lib/743d6dd4-ea73-4d4b-bbc7-d6de2ba45b73/scratchpad/mtk11-env` (temporary). `adapters.jl` `include`s the prototype's `BGProto.jl` by absolute path, from the `prototype/connection-mechanism` worktree.

```julia
# Adapters between BGProto's `Power` port and MTKSL connectors.
# Convention: an adapter is an ordinary bond-graph ELEMENT. The MTKSL network hanging off
# its MTKSL connector *is* the element at the end of the bond: `power.f` is the flow into
# it, `power.e` the effort across it. MTKSL through-variables are "into the adapter", so they
# are negated; the across-variable maps to the effort (electrical, hydraulic) or, for
# mechanics (mobility analogy), its rate maps to the flow.
include("/Users/fctc/me/doutorado/dev/bg/bg_lib-prototype/prototype/connection_mechanism/BGProto.jl")
using .BGProto
using ModelingToolkit
import ModelingToolkitStandardLibrary as MSL
const El = MSL.Electrical
const Rot = MSL.Mechanical.Rotational
const Tr = MSL.Mechanical.Translational
const TP = MSL.Mechanical.TranslationalPosition
const Hy = MSL.Hydraulic.IsothermalCompressible
const Bl = MSL.Blocks

function ElectricalAdapter(; name)          # two pins, like an MTKSL OnePort
    @named power = Power(); @named p = El.Pin(); @named n = El.Pin()
    System([power.e ~ p.v - n.v, power.f ~ -p.i, 0 ~ p.i + n.i], t; name, systems = [power, p, n])
end
function TranslationalAdapter(; name)       # velocity-based MechanicalPort (v, f)
    @named power = Power(); @named flange = Tr.MechanicalPort()
    System([power.e ~ -flange.f, power.f ~ flange.v], t; name, systems = [power, flange])
end
function RotationalAdapter(; name)          # position-based Flange (phi, tau), grounded support
    @named power = Power(); @named flange = Rot.Flange()
    System([power.e ~ -flange.tau, D(flange.phi) ~ power.f], t; name, systems = [power, flange])
end
function TranslationalPositionAdapter(; name)   # position-based Flange (s, f)
    @named power = Power(); @named flange = TP.Flange()
    System([power.e ~ -flange.f, D(flange.s) ~ power.f], t; name, systems = [power, flange])
end
function HydraulicAdapter(; name)           # (p, dm); BG flow = volume flow at reference density
    @named power = Power(); @named port = Hy.HydraulicPort()
    System([power.e ~ port.p, power.f ~ -port.dm / port.ρ], t; name, systems = [power, port])
end

# Causal blocks: a source whose effort/flow is a RealInput, a GY whose ratio is a RealInput,
# and an effort/flow sensor that writes a RealOutput.
function SeIn(; name)
    @named power = Power(); @named u = Bl.RealInput()
    System([power.e ~ u.u], t; name, systems = [power, u])
end
function mGY(; name)
    @named p1 = Power(); @named p2 = Power(); @named g = Bl.RealInput()
    System([p1.e ~ -g.u * p2.f, p2.e ~ g.u * p1.f], t; name, systems = [p1, p2, g])
end
function FlowSensor(; name)                 # zero-effort element: sits on a 1-junction
    @named power = Power(); @named y = Bl.RealOutput()
    System([power.e ~ 0, y.u ~ power.f], t; name, systems = [power, y])
end

using OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock
const TS = range(0.0, 10.0, length = 201)
simulate(model, u0 = []; alg = Rodas5P(), abstol = 1e-10) = (sys = mtkcompile(model);
    (sys, solve(ODEProblem(sys, u0, (0.0, 10.0)), alg; abstol, reltol = 1e-10)))
maxdiff(s1, x1, s2, x2) = maximum(abs, s1(TS; idxs = x1) .- s2(TS; idxs = x2))
```

## Appendix B: experiment scripts

### `coupling/probe_direct.jl`

```julia
# Does a direct `connect(Power, <MTKSL connector>)` work, error, or silently mis-pair?
include(joinpath(@__DIR__, "adapters.jl"))
function try_connect(label, other)
    @named bg = Damper(c = 1.0)
    try
        @named m = System([connect(bg.power, other)], t; systems = [bg, ModelingToolkit.get_parent(other) === nothing ? other : other])
        mtkcompile(m)
        println(rpad(label, 34), "COMPILED (no error)")
    catch err
        println(rpad(label, 34), "ERROR: ", first(split(sprint(showerror, err), '\n')))
    end
end
try_connect("Electrical.Pin (v, i)", El.Pin(name = :pin))
try_connect("Translational.MechanicalPort (v,f)", Tr.MechanicalPort(name = :mp))
try_connect("Rotational.Flange (phi, tau)", Rot.Flange(name = :fl))
try_connect("TranslationalPosition.Flange (s,f)", TP.Flange(name = :fs))
try_connect("Hydraulic HydraulicPort (p, dm)", Hy.HydraulicPort(name = :hp))
try_connect("Blocks.RealInput (u)", Bl.RealInput(name = :ri))

# Same names, same order as Pin -> silent positional pairing? (a Power renamed to v, i)
@connector function VIPower(; name)
    vars = @variables v(t) [guess = 0.0] i(t) [connect = Flow, guess = 0.0]
    System(Equation[], t, vars, []; name)
end
@named r = El.Resistor(R = 2.0); @named src = El.Voltage(); @named c = Bl.Constant(k = 3.0)
@named gnd = El.Ground(); @named port = VIPower()
@named bgR = System([port.v ~ 5.0 * port.i], t; systems = [port])   # a grounded BG-style R=5
@named m = System([connect(c.output, src.V), connect(src.p, r.p), connect(r.n, bgR.port),
                   connect(src.n, gnd.g)], t; systems = [r, src, c, gnd, bgR])
sys = mtkcompile(m)
println("\nVIPower(v,i) wired to Pin: compiled; equations after expansion:")
foreach(println, equations(ModelingToolkit.expand_connections(m)))
sol = solve(ODEProblem(sys, [], (0.0, 1.0)), Rodas5P())
println("current through series R(2)+bgR(5) at 3V: ", sol[r.i][end], "  (expect 3/7 = ", 3/7, ")")
```

### `coupling/ex1_rc.jl`

```julia
# (1) Bond-graph RC fed by an MTKSL Voltage source (driven by Blocks.Sine), vs pure MTKSL RC.
include(joinpath(@__DIR__, "adapters.jl"))
const Rv, Cv = 2.0, 0.5

# Reference: pure MTKSL series RC.
@named sine = Bl.Sine(frequency = 0.3, amplitude = 3.0, offset = 1.0)
@named src = El.Voltage(); @named gnd = El.Ground()
@named res = El.Resistor(R = Rv); @named cap = El.Capacitor(C = Cv, v = 0.5)
@named ref = System([connect(sine.output, src.V), connect(src.p, res.p), connect(res.n, cap.p),
                     connect(cap.n, src.n, gnd.g)], t; systems = [sine, src, gnd, res, cap])
_, sref = simulate(ref)

# (1a) R and C as bond-graph elements on a 1-junction (series), the source stays MTKSL.
# The adapter takes the place an Se would have: Junction1(:src, ...) gives it sign +1.
let
    @named sine = Bl.Sine(frequency = 0.3, amplitude = 3.0, offset = 1.0)
    @named src = El.Voltage(); @named gnd = El.Ground(); @named ea = ElectricalAdapter()
    @named R = Damper(c = Rv); @named C = Spring(k = 1 / Cv, x = 0.5 * Cv)   # q0 = C v0
    @named rc = Junction1(:src, [-1, R], [-1, C])
    @named m = System([connect(sine.output, src.V), connect(src.p, ea.p), connect(src.n, ea.n, gnd.g),
                       connect(ea.power, rc.src)], t; systems = [sine, src, gnd, ea, rc])
    sys, s = simulate(m)
    println("(1a) BG series RC via ElectricalAdapter: unknowns = ", unknowns(sys))
    println("     max |Δ v_C| = ", maxdiff(s, m.rc.C.power.e, sref, ref.cap.v),
            "   max |Δ i| = ", maxdiff(s, m.rc.R.power.f, sref, ref.res.i))
end

# (1b) Mixed: MTKSL resistor in the circuit, only the capacitor is bond-graph (node to ground).
let
    @named sine = Bl.Sine(frequency = 0.3, amplitude = 3.0, offset = 1.0)
    @named src = El.Voltage(); @named gnd = El.Ground(); @named res = El.Resistor(R = Rv)
    @named ea = ElectricalAdapter(); @named C = Spring(k = 1 / Cv, x = 0.5 * Cv)
    @named m = System([connect(sine.output, src.V), connect(src.p, res.p), connect(res.n, ea.p),
                       connect(ea.n, src.n, gnd.g), connect(ea.power, C.power)], t;
                      systems = [sine, src, gnd, res, ea, C])
    sys, s = simulate(m)
    println("(1b) MTKSL R + BG C via adapter: unknowns = ", unknowns(sys))
    println("     max |Δ v_C| = ", maxdiff(s, m.C.power.e, sref, ref.cap.v),
            "   max |Δ i| = ", maxdiff(s, m.res.i, sref, ref.res.i))
end

# (1c) Same as (1a) but the source is a BG Se driven by the Sine block through a RealInput.
let
    @named sine = Bl.Sine(frequency = 0.3, amplitude = 3.0, offset = 1.0)
    @named V = SeIn()
    @named R = Damper(c = Rv); @named C = Spring(k = 1 / Cv, x = 0.5 * Cv)
    @named rc = Junction1(:src, [-1, R], [-1, C])
    @named m = System([connect(sine.output, V.u), connect(V.power, rc.src)], t; systems = [sine, V, rc])
    sys, s = simulate(m)
    println("(1c) BG RC + SeIn(RealInput) driven by Blocks.Sine: max |Δ v_C| = ",
            maxdiff(s, m.rc.C.power.e, sref, ref.cap.v))
end
```

### `coupling/ex2_dcmotor.jl`

```julia
# (2) Bond-graph DC motor (electrical side + GY + load torque) driving an MTKSL rotational
#     Inertia + Damper on the shaft, vs the pure bond-graph DC motor of the prototype.
include(joinpath(@__DIR__, "adapters.jl"))

function bg_motor_parts()
    @named L = Mass(m = 0.5); @named R = Damper(c = 1.0); @named V = Se(12.0); @named τ = Se(1.0)
    @named je = Junction1(V, [-1, R], [-1, L], [-1, :gy])
    @named gy = GY(g = 0.01)
    (; je, gy, τ)
end

# Reference: pure bond graph (prototype run.jl's dc_motor).
let (; je, gy, τ) = bg_motor_parts()
    @named J = Mass(m = 0.01); @named c = Damper(c = 0.1)
    @named jm = Junction1(τ, [-1, c], [-1, J], :gy)
    @named ref = System([connect(je.gy, gy.p1), connect(gy.p2, jm.gy)], t; systems = [je, gy, jm])
    global ref = ref
    global (_, sref) = simulate(ref)
end

# (2a) Mechanical load in MTKSL: shaft = RotationalAdapter -> Inertia -> Damper -> Fixed.
let (; je, gy, τ) = bg_motor_parts()
    @named jm = Junction1(τ, :gy, [-1, :shaft])       # the adapter goes where J and c were (sign -1)
    @named ra = RotationalAdapter()
    @named inertia = Rot.Inertia(J = 0.01, phi = 0.0, w = 0.0)
    @named damper = Rot.Damper(d = 0.1); @named fixed = Rot.Fixed()
    @named m = System([connect(je.gy, gy.p1), connect(gy.p2, jm.gy), connect(jm.shaft, ra.power),
                       connect(ra.flange, inertia.flange_a), connect(inertia.flange_b, damper.flange_b),
                       connect(damper.flange_a, fixed.flange)], t;
                      systems = [je, gy, jm, ra, inertia, damper, fixed])
    sys, s = simulate(m)
    println("(2a) BG motor -> MTKSL Inertia+Damper: unknowns = ", string.(unknowns(sys)))
    println("     max |Δ ω| = ", maxdiff(s, m.inertia.w, sref, ref.jm.J.power.f),
            "   max |Δ i| = ", maxdiff(s, m.je.L.power.f, sref, ref.je.L.power.f))
end

# (2b) The other way round: MTKSL electrical circuit (Voltage, Resistor, Inductor) feeding the
#      BG GY + BG mechanical side through an ElectricalAdapter.
let
    @named k = Bl.Constant(k = 12.0); @named src = El.Voltage(); @named gnd = El.Ground()
    @named res = El.Resistor(R = 1.0); @named ind = El.Inductor(L = 0.5, i = 0.0)
    @named ea = ElectricalAdapter(); @named gy = GY(g = 0.01)
    @named J = Mass(m = 0.01); @named c = Damper(c = 0.1); @named τ = Se(1.0)
    @named jm = Junction1(τ, [-1, c], [-1, J], :gy)
    @named m = System([connect(k.output, src.V), connect(src.p, res.p), connect(res.n, ind.p),
                       connect(ind.n, ea.p), connect(ea.n, src.n, gnd.g),
                       connect(ea.power, gy.p1), connect(gy.p2, jm.gy)], t;
                      systems = [k, src, gnd, res, ind, ea, gy, jm])
    sys, s = simulate(m)
    # GY.p1 is a two-port's port: its flow is into the GY, i.e. the adapter is element-like
    # from the GY's view; here power flows adapter -> GY, so check the signs against the reference.
    println("(2b) MTKSL circuit -> adapter -> BG GY/mech: max |Δ ω| = ",
            maxdiff(s, m.jm.J.power.f, sref, ref.jm.J.power.f),
            "   max |Δ i| = ", maxdiff(s, m.ind.i, sref, ref.je.L.power.f))
end

# (2c) Pure MTKSL DC motor (EMF) as an independent sign-convention check of the BG reference.
let
    @named k = Bl.Constant(k = 12.0); @named src = El.Voltage(); @named gnd = El.Ground()
    @named res = El.Resistor(R = 1.0); @named ind = El.Inductor(L = 0.5, i = 0.0)
    @named emf = El.EMF(k = 0.01); @named fixed = Rot.Fixed()
    @named inertia = Rot.Inertia(J = 0.01, phi = 0.0, w = 0.0); @named damper = Rot.Damper(d = 0.1)
    @named load = Rot.ConstantTorque(tau_constant = 1.0)
    @named m = System([connect(k.output, src.V), connect(src.p, res.p), connect(res.n, ind.p),
                       connect(ind.n, emf.p), connect(emf.n, src.n, gnd.g),
                       connect(emf.flange, inertia.flange_a, load.flange), connect(emf.support, fixed.flange),
                       connect(inertia.flange_b, damper.flange_b), connect(damper.flange_a, fixed.flange)], t;
                      systems = [k, src, gnd, res, ind, emf, fixed, inertia, damper, load])
    sys, s = simulate(m)
    println("(2c) pure MTKSL motor (EMF) vs pure BG motor: max |Δ ω| = ",
            maxdiff(s, m.inertia.w, sref, ref.jm.J.power.f),
            "   max |Δ i| = ", maxdiff(s, m.ind.i, sref, ref.je.L.power.f))
end
```

### `coupling/ex3_translational.jl`

```julia
# (3) Bond-graph mass-spring-damper (forced by a Blocks.Sine via SeIn) coupled through an
#     MTKSL spring to an MTKSL mass, vs a pure MTKSL Translational model.
include(joinpath(@__DIR__, "adapters.jl"))
const m1v, k1v, d1v, k2v, m2v, x0 = 1.0, 1.0, 0.3, 2.0, 0.5, 0.5
sine() = Bl.Sine(frequency = 0.2, amplitude = 1.0, name = :sine)

# Reference: pure MTKSL, velocity-based Translational.
@named m1 = Tr.Mass(m = m1v); @named k1 = Tr.Spring(k = k1v); @named d1 = Tr.Damper(d = d1v)
@named k2 = Tr.Spring(k = k2v); @named m2 = Tr.Mass(m = m2v)
@named fx = Tr.Fixed(); @named F = Tr.Force(); sn = sine()
@named ref = System([connect(sn.output, F.f), connect(F.flange, m1.flange, k1.flange_a, d1.flange_a, k2.flange_a),
                     connect(k1.flange_b, d1.flange_b, fx.flange), connect(k2.flange_b, m2.flange)], t;
                    systems = [m1, k1, d1, k2, m2, fx, F, sn])
_, sref = simulate(ref, [m1.v => 0, m1.s => 0, m2.v => 0, m2.s => 0, k1.delta_s => x0, k2.delta_s => 0])

function bg_side()   # BG 1-DoF: mass, spring (x0), damper, input force, free port to the outside
    @named m = Mass(m = m1v); @named s = Spring(k = k1v, x = x0); @named d = Damper(c = d1v)
    @named f = SeIn(); sn = sine()
    @named j = Junction1([-1, m], [-1, s], [-1, d], f, [-1, :ext])   # adapter: sign -1, like m/s/d
    (; j, sn)
end

# (3a) velocity-based MechanicalPort through TranslationalAdapter
let (; j, sn) = bg_side()
    @named ta = TranslationalAdapter(); @named k2 = Tr.Spring(k = k2v); @named m2 = Tr.Mass(m = m2v)
    @named mdl = System([connect(sn.output, j.f.u), connect(j.ext, ta.power), connect(ta.flange, k2.flange_a),
                         connect(k2.flange_b, m2.flange)], t; systems = [j, sn, ta, k2, m2])
    sys, s = simulate(mdl, [m2.v => 0, m2.s => 0, k2.delta_s => 0])
    println("(3a) BG msd -> TranslationalAdapter -> MTKSL spring+mass: unknowns = ", string.(unknowns(sys)))
    println("     max |Δ v1| = ", maxdiff(s, mdl.j.m.power.f, sref, ref.m1.v),
            "   max |Δ v2| = ", maxdiff(s, mdl.m2.v, sref, ref.m2.v),
            "   max |Δ x1| = ", maxdiff(s, mdl.j.s.q, sref, ref.k1.delta_s))
end

# (3b) position-based TranslationalPosition.Flange (s, f) through TranslationalPositionAdapter
let (; j, sn) = bg_side()
    @named tpa = TranslationalPositionAdapter(); @named k2 = TP.Spring(k = k2v); @named m2 = TP.Mass(m = m2v, s = 0.0, v = 0.0)
    @named mdl = System([connect(sn.output, j.f.u), connect(j.ext, tpa.power), connect(tpa.flange, k2.flange_a),
                         connect(k2.flange_b, m2.flange)], t; systems = [j, sn, tpa, k2, m2])
    sys, s = simulate(mdl, [tpa.flange.s => 0.0])
    println("(3b) BG msd -> TranslationalPositionAdapter -> TP spring+mass: unknowns = ", string.(unknowns(sys)))
    println("     max |Δ v1| = ", maxdiff(s, mdl.j.m.power.f, sref, ref.m1.v),
            "   max |Δ v2| = ", maxdiff(s, mdl.m2.v, sref, ref.m2.v))
end
```

### `coupling/ex4_mgy.jl`

```julia
# (4) GY ratio from MTKSL Blocks: RealInput-modulated GY (mGY) vs the expression form GY(g = expr).
include(joinpath(@__DIR__, "adapters.jl"))

# The prototype's DC motor with a pluggable gyrator `gy`, extra subsystems and equations.
function motor(gy, extra_sys = [], extra_eqs = Equation[]; sensor = false)
    @named L = Mass(m = 0.5); @named R = Damper(c = 1.0); @named V = Se(12.0); @named τ = Se(1.0)
    @named J = Mass(m = 0.01); @named c = Damper(c = 0.1)
    @named je = Junction1(V, [-1, R], [-1, L], [-1, :gy])
    if sensor   # FlowSensor = zero-effort element on the shaft 1-junction, outputs ω
        @named ω = FlowSensor()
        @named jm = Junction1(τ, [-1, c], [-1, J], [-1, ω], :gy)
    else
        @named jm = Junction1(τ, [-1, c], [-1, J], :gy)
    end
    System([connect(je.gy, gy.p1); connect(gy.p2, jm.gy); extra_eqs], t;
           systems = [je, gy, jm, extra_sys...], name = :motor)
end
ω(s, m) = s(TS; idxs = m.jm.J.power.f)
cmp(label, a, b) = println(rpad(label, 66), "max |Δ ω| = ", maximum(abs, a .- b))

m0 = motor(GY(g = 0.01, name = :gy)); _, s0 = simulate(m0)

# (4a) constant ratio: Blocks.Constant -> mGY.g, vs GY(g = 0.01)
let @named k = Bl.Constant(k = 0.01); @named gy = mGY()
    m = motor(gy, [k], [connect(k.output, gy.g)]); _, s = simulate(m)
    cmp("(4a) Blocks.Constant -> mGY.g   vs GY(g = 0.01)", ω(s, m), ω(s0, m0))
end

# (4b) time-varying ratio: Blocks.Sine -> mGY.g, vs GY(g = 0.01 + 0.005 sin(2π·0.5 t))
mr = motor(GY(g = 0.01 + 0.005 * sin(2π * 0.5 * t), name = :gy)); _, sr = simulate(mr)
cmp("(sanity) GY(g = expr of t)   vs unmodulated GY(0.01)", ω(sr, mr), ω(s0, m0))
let @named sn = Bl.Sine(frequency = 0.5, amplitude = 0.005, offset = 0.01); @named gy = mGY()
    m = motor(gy, [sn], [connect(sn.output, gy.g)]); _, s = simulate(m)
    cmp("(4b) Blocks.Sine -> mGY.g   vs GY(g = expr of t)", ω(s, m), ω(sr, mr))
end
# (4b') no new component: the ratio expression reads the block's output variable directly.
let @named sn = Bl.Sine(frequency = 0.5, amplitude = 0.005, offset = 0.01)
    @named gy = GY(g = sn.output.u)          # @named wraps symbolic kwargs in ParentScope
    m = motor(gy, [sn]); _, s = simulate(m)
    cmp("(4b') @named gy = GY(g = sn.output.u)   vs GY(g = expr of t)", ω(s, m), ω(sr, mr))
    try
        simulate(motor(GY(g = sn.output.u, name = :gy), [sn]))   # same, without @named
        println("(4b'') GY(g = sn.output.u, name = :gy): compiled")
    catch err
        println("(4b'') GY(g = sn.output.u, name = :gy) without @named: ", first(split(sprint(showerror, err), '\n')))
    end
end

# (4c) closed loop through causal blocks: ω -> FlowSensor -> Gain(0.001) -> Add(+0.01) -> mGY.g,
#      vs the prototype's mGY expression g = 0.01 + 0.001 ω (run.jl, modulated = true).
let @named gain = Bl.Gain(k = 0.001); @named off = Bl.Constant(k = 0.01); @named add = Bl.Add()
    @named gy = mGY()
    inner = motor(gy; sensor = true)
    @named top = System([connect(inner.jm.ω.y, gain.input), connect(gain.output, add.input1),
                         connect(off.output, add.input2), connect(add.output, inner.gy.g)], t;
                        systems = [inner, gain, off, add])
    sys, s = simulate(top)
    @named L = Mass(m = 0.5); @named R = Damper(c = 1.0); @named V = Se(12.0); @named τ = Se(1.0)
    @named J = Mass(m = 0.01); @named c = Damper(c = 0.1)
    @named je = Junction1(V, [-1, R], [-1, L], [-1, :gy]); @named jm = Junction1(τ, [-1, c], [-1, J], :gy)
    @named gyr = GY(g = 0.01 + 0.001 * jm.J.power.f)
    @named ref = System([connect(je.gy, gyr.p1), connect(gyr.p2, jm.gy)], t; systems = [je, gyr, jm])
    _, sref = simulate(ref)
    cmp("(4c) FlowSensor->Gain->Add->mGY.g   vs GY(g = 0.01+0.001ω)",
        s(TS; idxs = top.motor.jm.J.power.f), sref(TS; idxs = ref.jm.J.power.f))
    println("     unknowns: ", string.(unknowns(sys)))
end
```

### `coupling/ex5_hydraulic.jl`

```julia
# (5) Hydraulic: MTKSL pressure / mass-flow sources (with HydraulicFluid) feeding bond-graph
#     R and C through HydraulicAdapter, vs pure bond-graph equivalents (Se / Sf).
include(joinpath(@__DIR__, "adapters.jl"))
const AT = 1e-18   # states are volumes q ~ 1e-6 m³: abstol must sit far below them
const p0, Rh, Ch, ρ = 2.0e5, 1.0e9, 2.0e-9, 997.0   # Pa, Pa·s/m³, m³/Pa, kg/m³ (water default)

# (5a) MTKSL Pressure source driven by Blocks.Constant -> adapter -> BG series R + C
@named Rr = Damper(c = Rh); @named Cr = Spring(k = 1 / Ch)
@named Se0 = Se(p0); @named ref_a = Junction1(Se0, [-1, Rr], [-1, Cr])
_, sa_ref = simulate(ref_a; abstol = AT)
let
    @named k = Bl.Constant(k = p0); @named src = Hy.Pressure(); @named fluid = Hy.HydraulicFluid(density = ρ)
    @named ha = HydraulicAdapter(); @named R = Damper(c = Rh); @named C = Spring(k = 1 / Ch)
    @named rc = Junction1(:src, [-1, R], [-1, C])
    @named m = System([connect(k.output, src.p), connect(fluid, src.port), connect(src.port, ha.port), connect(ha.power, rc.src)], t;
                      systems = [k, src, fluid, ha, rc])
    sys, s = simulate(m; abstol = AT)
    println("(5a) MTKSL Pressure -> HydraulicAdapter -> BG R+C: unknowns = ", string.(unknowns(sys)))
    println("     max |Δ p_C| = ", maxdiff(s, m.rc.C.power.e, sa_ref, ref_a.Cr.power.e),
            "   max |Δ Q| = ", maxdiff(s, m.rc.R.power.f, sa_ref, ref_a.Rr.power.f),
            "   check dm = -ρQ: ", maximum(abs, s(TS; idxs = m.src.port.dm) .+ ρ .* s(TS; idxs = m.rc.R.power.f)))
end

# (5b) MTKSL MassFlow source (dm from a Blocks.Sine) -> adapter -> BG R ∥ C (0-junction), vs BG Sf(dm/ρ)
const dmA = 1.0e-3   # kg/s amplitude
@named R2 = Damper(c = Rh); @named C2 = Spring(k = 1 / Ch); @named Sf0 = Sf(dmA / ρ * sin(2π * 0.2 * t))
@named ref_b = Junction0(Sf0, R2, C2)
_, sb_ref = simulate(ref_b; abstol = AT)
let
    @named sn = Bl.Sine(frequency = 0.2, amplitude = dmA); @named src = Hy.MassFlow()
    @named fluid = Hy.HydraulicFluid(density = ρ); @named ha = HydraulicAdapter()
    @named R = Damper(c = Rh); @named C = Spring(k = 1 / Ch)
    @named m = System([connect(sn.output, src.dm), connect(fluid, src.port), connect(src.port, ha.port),
                       connect(ha.power, R.power, C.power)], t; systems = [sn, src, fluid, ha, R, C])
    sys, s = simulate(m; abstol = AT)
    println("(5b) MTKSL MassFlow -> HydraulicAdapter -> BG R∥C: max |Δ p| = ",
            maxdiff(s, m.C.power.e, sb_ref, ref_b.C2.power.e), "  (p scale ~ ", maximum(abs, sb_ref(TS; idxs = ref_b.C2.power.e)), ")")
end
```

### `coupling/probe_init.jl`

```julia
# Position-based adapters add an absolute angle/position state. What if nobody gives it a value?
include(joinpath(@__DIR__, "adapters.jl"))
@named τ = Se(1.0); @named jm = Junction1(τ, [-1, :shaft])
@named ra = RotationalAdapter()
@named inertia = Rot.Inertia(J = 0.01)              # no phi / w given
@named damper = Rot.Damper(d = 0.1); @named fixed = Rot.Fixed()
@named m = System([connect(jm.shaft, ra.power), connect(ra.flange, inertia.flange_a),
                   connect(inertia.flange_b, damper.flange_b), connect(damper.flange_a, fixed.flange)], t;
                  systems = [jm, ra, inertia, damper, fixed])
sys = mtkcompile(m); println("unknowns: ", string.(unknowns(sys)))
for u0 in ([], [inertia.w => 0.0], [inertia.w => 0.0, inertia.phi => 0.0])
    try
        s = solve(ODEProblem(sys, u0, (0.0, 10.0)), Rodas5P())
        println("u0 = ", u0, " -> ", s.retcode, ", ω(10) = ", s[inertia.w][end], " (expect 10(1-e^-100) ≈ 10)")
    catch err
        println("u0 = ", u0, " -> ERROR: ", first(split(sprint(showerror, err), '\n')))
    end
end
```

### `coupling/probe_tmodelica.jl`

```julia
# TranslationalPositionAdapter (TP.Flange) connected to a TranslationalModelica.Mass:
# Se(1) on a BG 1-junction pushing an MTKSL Modelica-style mass m = 2 from rest: v(t) = t/2.
include(joinpath(@__DIR__, "adapters.jl"))
const TM = MSL.Mechanical.TranslationalModelica
@named F = Se(1.0); @named j = Junction1(F, [-1, :ext]); @named a = TranslationalPositionAdapter()
@named mass = TM.Mass(m = 2.0, s = 0.0, v = 0.0)
@named m = System([connect(j.ext, a.power), connect(a.flange, mass.flange_a)], t; systems = [j, a, mass])
sys, s = simulate(m, [a.flange.s => 0.0])
println("TP adapter -> TranslationalModelica.Mass: max |v - t/2| = ", maximum(abs, s(TS; idxs = mass.v) .- TS ./ 2))
```

### `coupling/probe_scope2.jl`

```julia
include(joinpath(@__DIR__, "adapters.jl"))
function motor(gy, extra_sys = [])
    @named L = Mass(m = 0.5); @named R = Damper(c = 1.0); @named V = Se(12.0); @named τ = Se(1.0)
    @named J = Mass(m = 0.01); @named c = Damper(c = 0.1)
    @named je = Junction1(V, [-1, R], [-1, L], [-1, :gy]); @named jm = Junction1(τ, [-1, c], [-1, J], :gy)
    System([connect(je.gy, gy.p1); connect(gy.p2, jm.gy)], t; systems = [je, gy, jm, extra_sys...], name = :motor)
end
for (label, mk) in [
        "GY(g = sn.output.u, name = :gy)" => () -> (sn = Bl.Sine(frequency = 0.5, amplitude = 0.005, offset = 0.01, name = :sn); motor(GY(g = sn.output.u, name = :gy), [sn])),
        "@named gy = GY(g = sn.output.u)" => () -> (@named sn = Bl.Sine(frequency = 0.5, amplitude = 0.005, offset = 0.01); @named gy = GY(g = sn.output.u); motor(gy, [sn])),
        "name :gyx" => () -> (@named sn = Bl.Sine(frequency = 0.5, amplitude = 0.005, offset = 0.01); @named gyx = GY(g = sn.output.u); motor(gyx, [sn])),
    ]
    m = mk()
    eqs = equations(ModelingToolkit.expand_connections(m))
    println(label, " -> ", filter(e -> occursin("output", string(e)) && occursin("p1₊e", string(e.lhs)), eqs))
    try mtkcompile(m); println("   compiles") catch e; println("   FAILS: ", first(split(sprint(showerror, e), '\n'))) end
end
```

## Appendix C: output of the final run (`coupling/results.log`)

```text
=== probe_direct
Electrical.Pin (v, i)             ERROR: Different types of connectors are in one connection statement: <[:bg₊power, :pin]>
Translational.MechanicalPort (v,f)ERROR: Different types of connectors are in one connection statement: <[:bg₊power, :mp]>
Rotational.Flange (phi, tau)      ERROR: Different types of connectors are in one connection statement: <[:bg₊power, :fl]>
TranslationalPosition.Flange (s,f)ERROR: Different types of connectors are in one connection statement: <[:bg₊power, :fs]>
Hydraulic HydraulicPort (p, dm)   ERROR: Different types of connectors are in one connection statement: <[:bg₊power, :hp]>
Blocks.RealInput (u)              ERROR: Different types of connectors are in one connection statement: <[:bg₊power, :ri]>
current through series R(2)+bgR(5) at 3V: 0.42857142857142855  (expect 3/7 = 0.42857142857142855)
=== ex1_rc
(1a) BG series RC via ElectricalAdapter: unknowns = SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}[rc₊C₊q(t)]
     max |Δ v_C| = 3.2053604215320775e-10   max |Δ i| = 1.6026813209890634e-10
(1b) MTKSL R + BG C via adapter: unknowns = SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}[C₊q(t)]
     max |Δ v_C| = 3.4673730553436144e-10   max |Δ i| = 1.7336887481178564e-10
(1c) BG RC + SeIn(RealInput) driven by Blocks.Sine: max |Δ v_C| = 3.2053604215320775e-10
=== ex2_dcmotor
(2a) BG motor -> MTKSL Inertia+Damper: unknowns = ["je₊L₊p(t)", "inertia₊w(t)", "inertia₊phi(t)"]
     max |Δ ω| = 5.674145597822644e-10   max |Δ i| = 1.6436985106338398e-10
(2b) MTKSL circuit -> adapter -> BG GY/mech: max |Δ ω| = 2.937774468136922e-10   max |Δ i| = 1.467590493575699e-10
(2c) pure MTKSL motor (EMF) vs pure BG motor: max |Δ ω| = 5.748095333046876e-10   max |Δ i| = 1.3927348163633724e-10
=== ex3_translational
(3a) BG msd -> TranslationalAdapter -> MTKSL spring+mass: unknowns = ["k2₊flange_b₊v(t)", "m2₊s(t)", "k2₊delta_s(t)", "j₊s₊q(t)", "j₊m₊p(t)"]
     max |Δ v1| = 3.0540348028296194e-11   max |Δ v2| = 5.259503943477739e-11   max |Δ x1| = 1.0938583372421817e-11
(3b) BG msd -> TranslationalPositionAdapter -> TP spring+mass: unknowns = ["k2₊flange_b₊s(t)", "tpa₊flange₊s(t)", "m2₊v(t)", "j₊s₊q(t)", "j₊m₊p(t)"]
     max |Δ v1| = 3.800526560127082e-11   max |Δ v2| = 7.658756961959057e-11
=== ex4_mgy
(4a) Blocks.Constant -> mGY.g   vs GY(g = 0.01)                   max |Δ ω| = 0.0
(sanity) GY(g = expr of t)   vs unmodulated GY(0.01)              max |Δ ω| = 0.56630349505377
(4b) Blocks.Sine -> mGY.g   vs GY(g = expr of t)                  max |Δ ω| = 0.0
(4b') @named gy = GY(g = sn.output.u)   vs GY(g = expr of t)      max |Δ ω| = 0.0
(4b'') GY(g = sn.output.u, name = :gy) without @named: ExtraVariablesSystemException: The system is unbalanced. There are 13 highest order derivative variables and 12 equations.
(4c) FlowSensor->Gain->Add->mGY.g   vs GY(g = 0.01+0.001ω)        max |Δ ω| = 7.691625114603085e-13
     unknowns: ["motor₊jm₊J₊p(t)", "motor₊je₊L₊p(t)"]
=== ex5_hydraulic
(5a) MTKSL Pressure -> HydraulicAdapter -> BG R+C: unknowns = ["rc₊C₊q(t)"]
     max |Δ p_C| = 0.0   max |Δ Q| = 0.0   check dm = -ρQ: 2.7755575615628914e-17
(5b) MTKSL MassFlow -> HydraulicAdapter -> BG R∥C: max |Δ p| = 1.9586195776355453e-8  (p scale ~ 488.5317146527688)
=== probe_init
unknowns: ["inertia₊w(t)", "inertia₊phi(t)"]
u0 = Any[] -> ERROR: Initialization incomplete. Not all of the state variables of the
u0 = Pair{Num, Float64}[inertia₊w(t) => 0.0] -> ERROR: Initialization incomplete. Not all of the state variables of the
u0 = Pair{Num, Float64}[inertia₊w(t) => 0.0, inertia₊phi(t) => 0.0] -> Success, ω(10) = 10.000000000257105 (expect 10(1-e^-100) ≈ 10)
=== probe_tmodelica
TP adapter -> TranslationalModelica.Mass: max |v - t/2| = 3.907985046680551e-14
=== probe_scope2
   FAILS: ExtraVariablesSystemException: The system is unbalanced. There are 13 highest order derivative variables and 12 equations.
   compiles
   compiles
```
