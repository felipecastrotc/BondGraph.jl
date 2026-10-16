# Spec: BondGraphToolkit v0.3 on ModelingToolkit 11

Status: ready-for-agent

Sources: the map (`map.md`) and its resolved tickets 01–06 and 08–11 (`issues/`), [ADR 0001](../../docs/adr/0001-junctions-on-native-mtk-connections.md), `GLOSSARY.md`, and the literature notes (`docs/literature/README.md`). Ticket 07 (licence) is deferred and doesn't block this spec.

## Problem Statement

BondGraph.jl v0.2 is a bond-graph dialect of ModelingToolkit (MTK). It works only on MTK 8.71 with a hand-kept Manifest.
- It overrides MTK internals to turn bonds into equations. MTK 10 rewrote those internals, and Julia 1.10+ refuses method overwriting during precompilation, so the package can't be installed with fresh dependency resolution.
- CI is red, and the package can't be registered or cited.

Its API has also drifted from the bond-graph literature:
- Electrical models build a capacitor from `Spring` and an inductor from `Mass`.
- Modulated two-ports need a tuple-returning `@named` form or a `coneqs` keyword.
- The inertia's state is its flow, so a nonlinear inertia can't be written.
- Several exports have no definition.
- There's no way to connect a bond graph to MTK standard-library (MTKSL) models.
- There's no readable way to get the model's equations for a report.
- The graph plot depends on removed MTK internals.

The maintainer wants to keep v0.2's ergonomics: bond-graph models are ordinary MTK systems, and custom nonlinear elements and modulated two-ports stay easy to write. The library must also follow bond-graph conventions exactly as the textbooks state them (Karnopp, Margolis & Rosenberg 2012; Borutzky 2011).

## Solution

v0.3 is renamed **BondGraphToolkit**. It runs on MTK 11 and Julia ≥ 1.10 and uses only MTK's public API.
- Bonds are ordinary MTK `connect`s between `PowerPort`s.
- Junctions are signed components that implement the textbook junction laws.
- Elements and two-ports are plain MTK systems.

A modeller builds a model the way they draw it: elements, junctions with each bond's half-arrow written as ±1, and two-ports with port 1 in and port 2 out. They then simulate with standard MTK (`mtkcompile` → `ODEProblem` → `solve`).

The library also provides:
- elements named by their textbook roles, with mechanical shorthands;
- detectors and signal-driven sources and moduli;
- adapters to MTKSL networks;
- `@oneport`/`@twoport` macros for custom elements;
- readable, renameable equations with LaTeX export;
- a bond-graph drawing (Graphviz, with TikZ export).

The physics of v0.2's examples doesn't change.

## User Stories

### Building models

1. As a modeller, I want to `using BondGraphToolkit` and get MTK re-exported, so that one import is enough to build and simulate a model.
2. As a modeller, I want bond-graph models to be ordinary MTK systems, so that I can compose them with any MTK model and use every MTK function on them.
3. As a modeller, I want to use MTK's own `t_nounits`/`D_nounits` in my custom elements, so that my expressions are interchangeable with any MTK model.
4. As a modeller, I want `Resistance(R = …)`, `Compliance(C = …, q = …)` and `Inertia(I = …, p = …)`, so that elements and parameters read as in the textbooks in every physical domain.
5. As a mechanical modeller, I want `Mass(m, u)`, `Spring(k, x)` and `Damper(c)` as shorthands, so that I can write mechanical models with stiffness and damping as in v0.2.
6. As a mechanical modeller, I want `Mass(m, u)` to start at the velocity `u`, so that I give initial conditions in the variable I know, while the state stays the momentum.
7. As a modeller, I want an inertia's state to be its momentum `p`, so that I can write a nonlinear inertia `f = Φ(p)` as the books do.
8. As a modeller, I want a compliance's state to be its displacement `q`, so that energy variables are consistent across C and I.
9. As a modeller, I want `Se(E)` and `Sf(F)` to take a number or any expression of time and system variables, so that forcing terms are one-liners.
10. As a modeller, I want `Junction1(se, [-1, m], [-1, d])` to mean bonds from `se` into the junction and from the junction out to `m` and `d`, so that the ±1 I write is exactly the sign in the junction sum.
11. As a modeller, I want `Junction0` and `Junction1` to honour the signs in the same way, so that one convention holds everywhere.
12. As a modeller, I want a 0-junction to give one shared effort and a signed flow sum, and a 1-junction one shared flow and a signed effort sum, so that the equations are the textbook ones.
13. As a modeller, I want element laws written in the variables of the bond as I drew it, so that the tool never reinterprets my drawing. A resistance with its bond drawn the other way acts as a negative resistor.
14. As a modeller, I want to give a junction named free ports (`:out`, `[-1, :out]`), so that I can bond junctions to other junctions and to two-ports with `connect`.
15. As a modeller, I want a plain `connect` of free ports to act as a 0-junction whose directions come from the signs at the far ends, so that multi-junction models need no extra component.
16. As a modeller, I want `TF(m = …)` with `e₁ = m e₂, f₂ = m f₁` and `GY(r = …)` with `e₁ = r f₂, e₂ = r f₁`, port 1 in and port 2 out, so that two-ports follow Karnopp Eqs. 3.9 and 3.11.
17. As a modeller, I want to pass an expression as the modulus, so that a modulated TF or GY (MTF, MGY) is the same constructor.
18. As a modeller, I want two-port ports named `p1`/`p2` and one-port ports named `power`, so that I can `connect` them predictably.
19. As a modeller, I want the connector type named `PowerPort`, so that it matches the glossary and doesn't clash with MTKSL's `Blocks.Power`.
20. As a modeller, I want a clear docstring saying an expression modulus needs `@named`, so that I don't fall into the scoping trap of calling a constructor with `name = …`.
21. As a modeller, I want each library element's docstring to show its conventional drawing direction, so that I know which sign to write at a junction.

### Signals and detectors

22. As a modeller, I want `Se()`/`Sf()` without a value to expose `input::RealInput`, so that a block diagram can drive a source (MSe, MSf).
23. As a modeller, I want `TF()`/`GY()` without a modulus to expose `input::RealInput`, so that a controller can modulate a two-port.
24. As a modeller, I want detectors `De` and `Df` with `output::RealOutput`, so that I can measure an effort or a flow without adding power to the bond.
25. As a control engineer, I want to close a loop from a detector through MTKSL blocks into a modulus, so that I can model feedback without leaving MTK.

### Coupling with the MTK standard library

26. As a modeller, I want `ElectricalAdapter` with pins `p`/`n`, so that I can attach an MTKSL circuit to a bond.
27. As a modeller, I want `TranslationalAdapter` (velocity-based) and `TranslationalPositionAdapter` (position-based, with an optional initial `s`), so that I can attach MTKSL translational models.
28. As a modeller, I want `RotationalAdapter` with an optional initial `phi`, so that I can attach MTKSL rotational models such as a load inertia.
29. As a modeller, I want `HydraulicAdapter` taking volume flow at the reference density, so that I can attach MTKSL isothermal-compressible networks.
30. As a modeller, I want every adapter to be a one-port with its law in bond variables, conventionally drawn into the adapter, so that it signs like any other element.
31. As a modeller, I want each adapter's docstring to state its sign rule and any gotcha (e.g. connect the hydraulic fluid in its own statement), so that coupling works the first time.

### Custom elements

32. As a modeller, I want `@oneport function MyC(; name, …) … [laws] end` to add the `power` port and wrap my laws in a system, so that a nonlinear R, C or I takes a few lines.
33. As a modeller, I want `@twoport` to add `p1` (bond in) and `p2` (bond out) and give me the bond variables `e1, f1, e2, f2`, so that I write a custom two-port's law exactly as in the book, without knowing that MTK measures port 2's flow inward.
34. As a modeller, I want the library's own elements built with these macros, so that custom elements are first-class and behave identically.
35. As a modeller, I want to write a multiport by hand as a plain MTK system with `PowerPort`s, documented with an example, so that I'm not limited to the macros.

### Equations, names and LaTeX

36. As a modeller, I want `bg_equations(sys)` to return the model's equations in terms of its states only, so that I can read the dynamics without junction internals.
37. As a modeller, I want `bg_equations(sys; rename = Dict(…))` to replace variable names, including inside derivatives (e.g. `D(s.q)` → `ẋ`, `m.p` → `p`), so that `m₊power₊f`-style names don't reach my report.
38. As a modeller, I want `latexify` available from the package, so that I can paste renamed equations into a report or paper.
39. As a modeller, I want documented guidance on aliasing a variable in the model (`x ~ s.q`), so that I can also read a renamed variable from the solution.

### Drawing the bond graph

40. As a modeller, I want a function that extracts the bond graph from a model (elements, junctions, two-ports, and bonds with their half-arrows), so that I can inspect or process the structure.
41. As a modeller, I want to render the bond graph with Graphviz, with half-arrow heads, `0`/`1` junction nodes and `type:parameter` labels laid out left to right, so that the picture looks like a textbook bond graph.
42. As an author, I want to export the drawing as TikZ in the style of the docs' figures, so that I can put it in a thesis or paper.
43. ~~As a modeller who doesn't draw graphs, I want the rendering dependency to be optional (a package extension), so that the core stays light.~~ Reversed by v0.3-impl ticket 11: Graphviz_jll is a regular dependency (batteries included).
44a. As a modeller, I want `bondgraph(model)` to show the bond graph as text, `bgplot(model)` to draw it like `plot(sol)`, and `savebg(path, model)` to export it by file extension, so that the drawing API reads like Plots (v0.3-impl ticket 11).

### Simulation and diagnostics

44. As a modeller, I want to simulate with `mtkcompile` → `ODEProblem` → `solve` and nothing package-specific, so that everything I know from MTK applies.
45. As a modeller, I want models where a storage element is in derivative causality (e.g. two inertias on one 1-junction) to compile, so that I'm not forced to merge elements by hand.
46. As a modeller, I want the docs to explain that MTK chooses states and reduces the index, and that a causal conflict surfaces as MTK's "structurally singular" error, so that I can interpret what I see.

### Examples and docs

47. As a reader, I want every element, port, macro and adapter to have a runnable docs example, so that I can learn each piece by copying working code.
48. As a reader, I want the six v0.2 examples (RC, DC motor, 1-, 2- and N-DoF mass–spring–damper, Duffing) rewritten for v0.3 with role names, so that they teach the current API.
49. As a reader, I want the N-DoF example drawn in the textbook orientation (`0 → sd`), so that the example matches the books.
50. As a reader, I want the Duffing example to keep its linear and cubic C on one 1-junction, so that it shows a custom element next to a library one, exactly as drawn.

### Packaging

51. As a maintainer, I want the package, module, main file, docs and tests renamed to `BondGraphToolkit` with the UUID kept, so that it can be registered.
52. As a maintainer, I want `ModelingToolkit = "11.45"` (caret), Julia ≥ 1.10, and direct dependencies only for what the code imports, so that fresh resolution works and bounds are honest.
53. As a maintainer, I want CI on Julia 1.10 and the latest 1.x to pass with no Manifest, so that the badge is green and registration can follow.
54. As a maintainer, I want dead code gone (`src/mtk.jl` overrides, bond metadata types, `simplifysys`, `substitute_dict`, `renamevars`, `addsubsys`, `Dq`, `isindependent`, `It`/DomainSets, `src/mbg.jl`, `generate_graph`), so that what remains is maintainable.

## Implementation Decisions

### Package

- Rename to `BondGraphToolkit` everywhere in code, docs and tests, keeping UUID `41e83443-…` and no alias (ticket 08). The GitHub repo rename happens at the release, after asking.
- Compat: `ModelingToolkit = "11.45"`, `julia = "1.10"`. Direct dependencies only for what the code imports, each with a `[compat]` bound:
  - ModelingToolkitStandardLibrary is a regular dependency (MIT, about 0.2 s load). It's needed for `RealInput`/`RealOutput` and the adapters.
  - Latexify is a direct dependency, because `latexify` is re-exported. MTK already loads it.
- Graph rendering uses Graphviz_jll as a regular dependency (v0.3-impl ticket 11; was a package extension).
- Licence: unchanged until ticket 07 is decided.
- CI: Julia `1.10` and `1`, no Manifest, no downgrade or nightly jobs.

### Sign convention (ticket 06, ADR 0001)

- The half-arrow is the positive direction of power. Both ends of a bond share one `e` and one `f`.
- Junctions: the shared variable is unsigned; the sum is signed, + for a bond into the junction and − for one out. This is the same at 0- and 1-junctions (Karnopp Eqs. 3.16–3.21).
- Element laws are written in the variables of their bond as drawn, with no built-in reorientation (a design choice; the books are silent on reversed one-ports).
- TF and GY are through-oriented: port 1 in, port 2 out.

### Junctions

- One signed junction component serves both kinds. The user-facing `Junction0(args…)`/`Junction1(args…)` build the component plus its members.
- Member syntax is unchanged from v0.2:
  - `x`: the bond points into the junction, +1.
  - `[-1, x]`: the bond points out of the junction.
  - `:name` or `[-1, :name]`: a free port, exposed as the junction system's own `PowerPort`.
- The junction converts MTK's port flow (positive into the component) to the bond flow. Prototype snippet, trimmed:

  ```julia
  # x = σ·c, where c is MTK's flow into the core port:
  #   element member: σ = -1 (the element's own variable is the bond flow)
  #   free port:      σ = s  (its sign)
  x = [σ * p.f for (σ, p) in zip(σs, ps)]
  kind == 0 ? [[p.e ~ ps[1].e for p in ps[2:end]]; 0 ~ sum(signs .* x)] :
              [[xk ~ x[1] for xk in x[2:end]]; 0 ~ sum(s * p.e for (s, p) in zip(signs, ps))]
  ```

- A junction records its members' signs and kinds as system metadata, so that the bond-graph extractor needs no MTK internals.

### Elements

- `PowerPort`: variables `e` and `f`, with `f` marked as MTK `Flow`.
- `Se(E)`, `Sf(F)`: the value is a number or an expression. Called with no value, they expose `input::RealInput`.
- `Resistance(R)`: `e = R f`.
- `Compliance(C, q)`: state `q`, with `e = q/C` and `q̇ = f`.
- `Inertia(I, p)`: state `p`, with `ṗ = e` and `f = p/I`.
- Shorthands, with v0.2 keywords and defaults:
  - `Mass(m, u)` is `Inertia(I = m, p = m·u)`.
  - `Spring(k, x)` is `Compliance(C = 1/k, q = x)`.
  - `Damper(c)` is `Resistance(R = c)`.

  They're exported and documented as mechanical shorthands. All examples use the role names.
- `De` and `Df`: zero flow or zero effort respectively on their bond, with `output.u` equal to the measured effort or flow.
- `TF(m)`: `e1 = m e2`, `f2 = m f1`. `GY(r)`: `e1 = r f2`, `e2 = r f1`. Both use bond variables. A modulus that is an expression gives the modulated form; no modulus exposes `input::RealInput`. There are no separate `MTF`/`MGY` names.
- Library elements don't apply `ParentScope` themselves. An expression modulus requires `@named`, and the docstrings say so.

### Macros

- `@oneport` adds `power::PowerPort` and wraps the returned laws in an MTK `System`. MTK 11 discovers unknowns and parameters itself, so v0.2's variable-sorting logic goes away.
- `@twoport` adds `p1`/`p2` and binds the bond variables `e1 = p1.e`, `f1 = p1.f`, `e2 = p2.e`, `f2 = −p2.f` before the user's laws.
- Both macros produce MTK components (`@component` semantics), and the library's one-ports and two-ports are written with them.

### Adapters (ticket 11)

The table below comes from research 11:

| Adapter | MTKSL connector(s) | Laws in bond variables |
|---|---|---|
| `ElectricalAdapter` | `Pin` `p`, `n` | `e = p.v − n.v`, `f = −p.i`, `p.i + n.i = 0` |
| `TranslationalAdapter` | `MechanicalPort` `flange` | `e = −flange.f`, `f = flange.v` |
| `TranslationalPositionAdapter` | `TranslationalPosition.Flange` `flange` | `e = −flange.f`, `D(flange.s) = f` |
| `RotationalAdapter` | `Rotational.Flange` `flange` | `e = −flange.tau`, `D(flange.phi) = f` |
| `HydraulicAdapter` | `HydraulicPort` `port` | `e = port.p`, `f = −port.dm/port.ρ` |

- Each adapter is conventionally drawn into the adapter. The table's signs were verified under the withdrawn convention A, so they must be re-verified against the MTKSL reference models under this spec's convention before release.
- The position-based adapters take `s = nothing` / `phi = nothing` for the initial value, with no hard default.
- A moving support (`use_support`) is out of scope.

### Equations and LaTeX

- `bg_equations(sys; rename = Dict())` takes a model (compiled or not) and returns MTK's full equations in terms of the states. It applies `rename` everywhere, including inside `Differential` (MTK's `substitute` doesn't reach there). The model is never modified.
- `latexify` is re-exported from Latexify.
- Measured on MTK 11.45: `latexify(full_equations(…))` is readable, while `latexify(compiled_system)` shows junction internals. The docs use the former, through `bg_equations`.

### Bond-graph drawing

- An exported extractor returns the bond graph as plain data:
  - nodes: each element with its type and parameter label, each junction with its kind, each two-port;
  - bonds: each with its from/to nodes per its half-arrow.
- The DOT writer draws bonds with half-arrow heads, junctions as `0`/`1`, labels as `type:name`, and lays out left to right in the textbook layout: the power path straight, each element beside its node.
- The TikZ writer uses Graphviz's node positions and emits a `.tex` file in the docs' figure style.
- Public API (v0.3-impl ticket 11): `bondgraph(sys)` is data shown as text; `bgplot(x; format = :svg | :tikz | :dot)` is the drawing or its source; `savebg(path, x)` exports `.svg .png .pdf .tex .dot`. No `plot`/`print` methods on MTK's `System` (type piracy).
- No causal strokes in v0.3.

### Removed

These go away:
- `src/mtk.jl` and every override of MTK internals;
- the bond metadata types and connection utilities;
- `Power`, `mGY`/`mTF` and the tuple `@named` form, `coneqs`;
- `simplifysys`, `reducedobs`, `substitute_dict`, the `structural_simplify` export;
- `renamevars`, `addsubsys`, `Dq`, `isindependent`, `It` and DomainSets;
- `generate_graph` and Plots/GraphRecipes;
- `src/mbg.jl`.

## Testing Decisions

- **Good tests check external behaviour only.** A test builds a model through the public API, then compiles, simulates or presents it, and checks the result against data, a textbook law or a reference model. Tests never inspect junction internals or the shape of the generated equations.
- **Seam 1: public API → `mtkcompile` → `solve`.**
  - Golden regression: the existing `test/regression` harness and CSVs. The six examples are rewritten for v0.3 and must reproduce the v0.2 data (`rtol 1e-6`). Only state names are remapped in the helpers, plus the sign of the N-DoF coupling springs' `q`.
  - Physics checks on small models:
    - junction power balance `Σ bᵢ eᵢ fᵢ = 0` along a trajectory, at both junction kinds;
    - a reversed `Resistance` delivers power (`q = e^{+t}` in the RC case);
    - `e1 f1 = e2 f2` for constant and modulated TF and GY;
    - `d(p²/2I)/dt = e·f` for `Inertia`;
    - two inertias on one 1-junction compile;
    - free-port junction-to-junction models give the textbook signs;
    - the `input`/`output` forms and `De`/`Df` match the value forms;
    - `@oneport`/`@twoport` custom elements match the built-ins;
    - the shorthands match the role-named elements.
  - Adapters: each matches a pure-MTKSL reference model (research 11's RC, DC motor, translational, rotational and hydraulic cases, plus a closed loop through `Df`).
- **Seam 2: presentation functions on a model.**
  - `bg_equations` with `rename`, including inside derivatives.
  - `latexify` of its output.
  - The extractor on a known model: nodes, kinds, labels, bond directions.
  - The DOT and TikZ text: half-arrow attributes, junction labels, one edge per bond.
  - Graphviz rendering: a smoke test in the extension only. No pixel comparisons.
- The unit tests of removed internals are deleted (`test_types`, `test_connection_utils`, `test_graph`, `test_bg_utils`, and the internal parts of `test_ports`/`test_system_utils`/`test_macros`).
- Prior art:
  - the golden harness (`test/regression/helpers.jl`, `generate_references.jl`);
  - the prototype's check script on branch `prototype/connection-mechanism` (12 checks, including the reversed R, the textbook N-DoF and mGY power balance);
  - research 11's coupling scripts.
- Docs examples run in the docs build. Every element, port, macro and adapter has one.

## Out of Scope

- A causality API. The SCAP diagnostic is a future option and never replaces MTK's equations by default.
- Causal strokes on the drawing.
- Multiport fields and a `@multiport` macro.
- Thermal and magnetic adapters. A moving support on mechanical adapters.
- Activity/power diagnostics, De/Df dualisation, sensitivity (incremental bond graphs).
- The licence decision and the MTK vs ModelingToolkitBase dependency (ticket 07, deferred).
- The documentation overhaul (roadmap step 4), beyond the per-feature examples and rewritten examples this spec requires.
- Registration in General and the software paper.
- 2D/multibond experiments.

## Further Notes

- Branch: `feature/mtk-upgrade` off `develop`, released as v0.3.0. At release, ask before renaming the GitHub repo.
- Follow-ups recorded elsewhere:
  - fix the DC motor and 2-DoF TikZ drawings to match the books (roadmap step 4);
  - fold `TODO.md` into tickets and delete it.
- Facts measured on MTK 11.45 during grilling:
  - two inertias on one 1-junction compile to one momentum state plus a dummy-derivative algebraic variable;
  - two `Sf` on one 1-junction raise "structurally singular", naming a port variable.

  Docs should show both.
- Literature references for every law above: `docs/literature/README.md`. Keep book statements separate from design choices in docstrings.
