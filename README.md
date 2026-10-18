<h1 align="center" margin=0px>
  <picture>
    <source srcset="./docs/src/assets/logo-h-dark.svg" media="(prefers-color-scheme: dark)">
    <img src="./docs/src/assets/logo-h.svg" width=70% alt="BondGraphToolkit.jl">
  </picture>
</h1>

<p align="center">
  <a href="https://www.felipecastrotc.com/BondGraph.jl/stable/"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Documentation (stable)"></a>
  <a href="https://www.felipecastrotc.com/BondGraph.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Documentation (dev)"></a>
  <a href="https://github.com/felipecastrotc/BondGraph.jl/actions/workflows/CI.yml?query=branch%3Amain"><img src="https://github.com/felipecastrotc/BondGraph.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="CI"></a>
  <a href="https://github.com/felipecastrotc/BondGraph.jl/releases/latest"><img src="https://img.shields.io/github/v/release/felipecastrotc/BondGraph.jl" alt="Latest release"></a>
  <img src="https://img.shields.io/badge/julia-%E2%89%A5%201.10-9558b2.svg" alt="Julia ≥ 1.10">
  <a href="LICENSE.md"><img src="https://img.shields.io/badge/license-Apache%202.0-blue.svg" alt="License: Apache 2.0"></a>
</p>

*Bond graphs in Julia, as ordinary ModelingToolkit systems.*

BondGraphToolkit.jl models dynamic systems with the bond-graph method, across mechanical, electrical, hydraulic and chemical domains. You build a model the way you draw it: elements, 0- and 1-junctions, transformers and gyrators, joined by bonds with a direction. The result is a plain [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl) system, so you can compose it with any MTK model and simulate it with the SciML solvers. The conventions follow the textbooks (Karnopp, Margolis & Rosenberg; Borutzky).

## Installation

The package is not yet in the General registry; add it from this repository:

```julia
pkg> add https://github.com/felipecastrotc/BondGraph.jl
pkg> add OrdinaryDiffEqTsit5   # an ODE solver for the examples
```

It requires Julia 1.10 or later and ModelingToolkit 11.

## Example

A capacitor discharging through a resistor: both share one voltage, so they hang off a 0-junction, with each bond pointing into the element.

```julia
using BondGraphToolkit
using OrdinaryDiffEqTsit5

@named C = Compliance(C = 1.0, q = 1.0)   # capacitance, initial charge
@named R = Resistance(R = 1.0)
@named model = Junction0(R, C)            # 0 ⇀ R, 0 ⇀ C

sys = mtkcompile(model)
sol = solve(ODEProblem(sys, [], (0.0, 10.0)), Tsit5())
sol(1.0; idxs = model.C.q)                # q₀ e^(-t/RC) = e⁻¹ ≈ 0.368

bgequations(model)                        # the model's equations, readable
bgplot(model)                             # its bond graph, drawn
```

## Features

- **Bonds with a direction.** Every element has an incoming port `power` and an outgoing one `pout`, so `connect(se.pout, tf.p1)` is `Se ⇀ TF`.
- **Textbook elements and junctions.** `Se`, `Sf`, `Resistance`, `Compliance`, `Inertia`, the detectors `De` and `Df`, signed `Junction0` and `Junction1`, and `TF` and `GY`, modulated by any expression or signal. Mechanical shorthands `Mass`, `Spring` and `Damper`.
- **Custom elements** in a few lines with `@oneport` and `@twoport`.
- **Signals and the MTK Standard Library.** Blocks drive sources and moduli, and adapters couple a bond graph to electrical, mechanical and hydraulic MTK networks.
- **Readable equations** with names of your choice, and LaTeX output (`bgequations`, `latexify`).
- **Drawings** in the textbook layout, as SVG, TikZ or DOT (`bgplot`, `bgsave`).

## Documentation

The [documentation](https://www.felipecastrotc.com/BondGraph.jl/stable/) covers the conventions, the API, and worked examples from mechanical systems to a DC-motor pendulum, an arterial Windkessel and a chemical reaction network, each checked against hand-derived equations. The changes in each version are in [NEWS.md](NEWS.md).

## Citing

If you use BondGraphToolkit.jl in your work, please cite it. GitHub's *Cite this repository* button gives the reference from [CITATION.cff](CITATION.cff).

## Contributing

Bug reports, questions and pull requests are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

Apache License 2.0; see [LICENSE.md](LICENSE.md).
