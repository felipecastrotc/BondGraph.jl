# BondGraphToolkit.jl release notes

## v0.3.0 (2026-10-18)

BondGraphToolkit 0.3.0 rebuilds the package on **ModelingToolkit 11** and **Julia ≥ 1.10**. It uses only MTK's public API, so it installs with a fresh dependency resolution: no pinned Manifest and no overridden MTK internals. A bond graph is an ordinary MTK system. You build it the way you draw it, then simulate it with `mtkcompile`, `ODEProblem` and `solve`.

This is a breaking release. The package was `BondGraph` (v0.2, MTK 8.71) and is now `BondGraphToolkit`. See *Upgrading from 0.2* below.

### Highlights

- **Bonds are MTK `connect`s with a direction.** Every element has an incoming port `power` and an outgoing port `pout`:
  - `connect(se.pout, tf.p1)` is `Se ⇀ TF`;
  - `connect(tf.p2, mass.power)` is `TF ⇀ I`;
  - a plain `connect` of several ports is a correctly signed 0-junction.
- **Textbook junctions.** `Junction0` and `Junction1` give one shared effort or flow and a signed sum (Karnopp, Margolis & Rosenberg). Members say which way each bond points:

  ```julia
  @named model = Junction1(se => :pout, R, C)   # Se ⇀ 1 ⇀ R, C
  ```

- **Elements named by their roles**: `Se`, `Sf`, `Resistance`, `Compliance`, `Inertia` (state: momentum `p`), the mechanical shorthands `Mass`, `Spring`, `Damper`, and the detectors `De`, `Df`.
- **Two-ports** `TF` and `GY`, modulated by any expression or signal.
- **Signals.** Sources and moduli read MTK Standard Library blocks, and detectors write to them.
- **Coupling with the MTK Standard Library**: electrical, translational, translational-position, rotational and hydraulic adapters. They conserve power whichever side supplies it.
- **Custom elements** with `@oneport` and `@twoport`.
- **Readable equations**: `bgequations(model; rename = …)`, with LaTeX output through `latexify`.
- **Drawings**: `BondGraph(model)` gives the bond graph as data, `bgplot(model)` draws it in the textbook layout (SVG, TikZ, DOT), and `bgsave("model.pdf", model)` saves it.
- **Documentation**: the conventions page, 19 worked examples (mechanical, electrical, electromechanical, control, robotics, biomedical, hydraulic, chemical), each checked against hand-derived equations or closed forms, and the v0.2 results kept as regression tests.

### Upgrading from 0.2

| 0.2 (`BondGraph`) | 0.3 (`BondGraphToolkit`) |
|---|---|
| `using BondGraph` | `using BondGraphToolkit` |
| `Junction1(se, [-1, m], [-1, d])` | `Junction1(se => :pout, m, d)`: a bare element is a bond into it; `x => :pout` is a bond out of it |
| free ports `[-1, :a]` / `:a` | `:a` (leaving the junction) / `:a => :in` (entering it) |
| `Spring` as a capacitor, `Mass` as an inductor | `Compliance`, `Inertia` (`Mass`, `Spring`, `Damper` remain as mechanical shorthands) |
| `Mass` state: velocity | `Inertia` state: momentum `p`; `Mass(m = …, v = …)` takes the initial velocity |
| `mTF`, `mGY`, `coneqs` | `TF(m = …)`, `GY(r = …)`: the modulus can be any expression or signal |
| `generate_graph` | `BondGraph(model)`, `bgplot(model)`, `bgsave(path, model)` |
| `structural_simplify` | `mtkcompile` (MTK 11) |
| `renamevars`, `simplifysys` | `bgequations(model; rename = Dict(…))` |

The physics of the 0.2 examples is unchanged: their results are kept as regression tests and still match.

### Requirements

Julia 1.10 or later, ModelingToolkit 11. Tested on Julia 1.10 and the latest release.

### Known limitations

- The library assigns no causal strokes; MTK picks the states and reduces the index. A state choice can be singular for some models; the docs show how to steer it with `state_priority`.
- A subsystem written as a component with several junctions is drawn as one node.
- Multiport fields (C-, I-, R-fields) aren't provided yet.


## v0.2.0 (2023-10-06)

`BondGraph` on ModelingToolkit 8.71: elements, signed junctions, modulated two-ports (`mTF`, `mGY`), `@oneport` and graph generation. The `vMTK-v8.71.2` tags are hotfixes pinned to that MTK version.
