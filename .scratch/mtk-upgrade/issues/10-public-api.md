# Decide the public API and its breaking changes

Type: grilling
Status: resolved

## Already decided (ticket 04)

- Junction-to-junction bonds via named free ports + explicit `connect`; `Mass` has momentum state `p`; `GY`/`TF` components replace `mGY(...; coneqs)` and the tuple `@named` form.

## Inputs

- [research/10-api-usage-inventory.md](../research/10-api-usage-inventory.md) (Codex): every export with signature, call sites in examples/tests/docs, usage patterns, mechanical names used for electrical elements (RC uses `Spring`/`Damper`, DC motor `Mass` for the inductor), dead exports (`renamevars`, `addsubsys`, `latexify`), and which symbols depend on removed MTK internals.

- [research/11-coupling-with-mtk-stdlib.md](../research/11-coupling-with-mtk-stdlib.md): adapters, `RealInput` forms, sensors, MTKSL dependency, and the `Power`/`Mass`/`Spring`/`Damper` name clashes with MTKSL.

## Question

What is the v0.3 public API: element names (`Mass`/`Spring`/`Damper` vs I/C/R per `GLOSSARY.md`), the `@named sys, con = mGY(...)` tuple form, the future of `@oneport` for custom elements, whether `generate_graph` (Plots/GraphRecipes) becomes a package extension, and the exports with no local definition (`renamevars`, `addsubsys`, `latexify`)?

## Answers so far (grilling, 2026-10-15)

- Q1 Element names: role names `Resistance`, `Compliance`, `Inertia` are the standard and are used by the library and examples. `Mass`/`Spring`/`Damper` stay as the maintainer's mechanical shorthands.
- Q2 Role-named elements take the book's parameters and states: `Resistance(R)`, `Compliance(C, q)`, `Inertia(I, p)`. The shorthands keep their v0.2 keywords (`Mass(m, u)`, `Spring(k, x)`, `Damper(c)`).
- Q3/Q14 Half-arrow at a junction: keep `[-1, x]` (bond out of the junction) and bare `x` (into); the number is the book's bᵢ in the junction sum. `-x` rejected: type piracy on MTK's `System`, and the sign belongs to the bond, not the element.
- Q4 Connector `PowerPort`; one-port port `power`; two-port ports `p1` (bond in), `p2` (bond out).
- Q5 One `TF(m = …)` and one `GY(r = …)` (Karnopp's modulus symbols); an expression gives the modulated form; expression moduli need `@named` (documented, no internal `ParentScope`).
- Q6 `Se()`/`Sf()`/`TF()`/`GY()` without a value expose `input::RealInput`; detectors `De`/`Df` with `output::RealOutput`; MTKSL in `[deps]`.
- Q7 Five MTKSL adapters (`Electrical`, `Translational`, `TranslationalPosition`, `Rotational`, `Hydraulic`), MTKSL port names, conventionally drawn into the adapter; signs re-verified against MTKSL references during implementation.
- Q8/Q15 Keep `@oneport`; add `@twoport` (`p1` in, `p2` out, bond variables `e1, f1, e2, f2`, so the law is written as in the book); `TF`/`GY` use it. Every element, port and macro gets a docs example (maintainer: a must).
- Q9/Q16 Drop `simplifysys` and `substitute_dict`; `full_equations` covers the substitution. Renaming variables for display (e.g. `m.p` → `p`, `D(s.q)` → `ẋ`) and LaTeX export of the equations are an ergonomics must (design in round 3).
- Q10 Drop `renamevars`, `addsubsys`, `latexify` (as dead exports), `Dq`, `isindependent`, `It`/DomainSets.
- Q11/Q17 Rebuild graph plotting: a public bond-graph extractor, Graphviz rendering (half-arrows, `0`/`1`, `type:parameter` labels), plus a TikZ export.
- Q12/Q18 Causality is left to MTK. TODO (future): SCAP as an optional diagnostic/alternative, never replacing MTK's equations by default.
- Q13 Duffing keeps its two C elements.
- Q19 `Mass`/`Spring`/`Damper` are exported with short docstrings as mechanical shorthands; all examples use the role names; `Mass(m, u)` starts at p₀ = m·u; `Damper` drops `u`.
- Q21 `bg_equations(sys; rename = Dict(...))` returns the full equations in terms of states with display names applied, also inside derivatives; the model is untouched. `latexify` is re-exported from Latexify. `sol[x]` access goes through an alias equation in the model (documented). Ships in v0.3 with a docs example.
- Q22 TikZ export ships in v0.3, using Graphviz's layout, in the docs' figure style.

Facts found (MTK 11.45, prototype): two I on one 1-junction compile to one momentum state plus a dummy-derivative algebraic; two Sf on one 1-junction raise "structurally singular"; `substitute` renames everywhere except inside `Differential`; `latexify(mtkcompile(m))` shows junction internals (`core.p1.f`), `latexify(full_equations(…))` is readable.

## Answer

The v0.3 public API (module `BondGraphToolkit`, `@reexport using ModelingToolkit`, time `t_nounits`/`D_nounits` from MTK):

- **Connector:** `PowerPort` (`e`, `f` with `Flow`).
- **One-ports:** `Se`, `Sf`, `Resistance(R)`, `Compliance(C, q)`, `Inertia(I, p)`; shorthands `Mass(m, u)`, `Spring(k, x)`, `Damper(c)`. Without a value, `Se()`/`Sf()` expose `input::RealInput`.
- **Detectors:** `De`, `Df` with `output::RealOutput`.
- **Two-ports:** `TF(m)`, `GY(r)`, port 1 in, port 2 out; an expression modulus gives MTF/MGY; no value exposes `input`.
- **Junctions:** `Junction0`, `Junction1` with members `x` (bond in, +1), `[-1, x]` (bond out), free ports `:name`/`[-1, :name]`.
- **Custom elements:** `@oneport`, `@twoport` (bond variables `e1, f1, e2, f2`).
- **MTKSL adapters:** `ElectricalAdapter`, `TranslationalAdapter`, `TranslationalPositionAdapter`, `RotationalAdapter`, `HydraulicAdapter`.
- **Presentation:** `bg_equations(sys; rename)`, `latexify`; a bond-graph extractor with Graphviz rendering and TikZ export.
- **Dropped:** `Power`, `mGY`/`mTF` and their tuple form, `simplifysys`, `substitute_dict`, `structural_simplify` export, `renamevars`, `addsubsys`, `Dq`, `isindependent`, `It`/DomainSets, `generate_graph` (replaced).

Every element, port, macro and adapter gets a docs example. Causality is left to MTK; SCAP as an optional diagnostic is a future TODO.
