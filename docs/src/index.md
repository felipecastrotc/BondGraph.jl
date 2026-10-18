# Welcome to Bond Graph Toolkit Documentation

*A bond graph toolkit for Julia.*

The `BondGraphToolkit.jl` is a Julia library designed to simplify the modeling and simulation of bond graphs systems. It offers a user-friendly syntax where engineers, researchers, and enthusiasts can express and analyze complex dynamic systems effortlessly. Built on top of the `ModelingToolkit.jl` ecosystem, the `BondGraphToolkit.jl` seamlessly integrates with existing tools, making it easy to use and customize.

## Package Features

The key features of the package include:

1. **Bond graphs as you draw them**: elements, 0- and 1-junctions and two-ports follow the textbook conventions (Karnopp, Margolis & Rosenberg), with each bond's half-arrow written as a `±1` sign.
2. **The standard elements**: sources (`Se`, `Sf`), resistance, compliance, inertia, transformers and gyrators (constant or modulated), and effort and flow detectors.
3. **Native ModelingToolkit models**: every element and model is an MTK system, simulated with `mtkcompile`, `ODEProblem` and `solve`, and coupled with any MTK model, including the MTK Standard Library through adapters and signal blocks.
4. **Custom elements**: nonlinear one-ports and two-ports in a few lines with `@oneport` and `@twoport`.
5. **Readable output**: the model's equations with your own variable names and in LaTeX, and drawings of the bond graph (SVG, or TikZ for a paper).

## Use Cases

- Control systems design and analysis: [DC motor speed control](@ref control-dc-motor-pi).
- Robotics and mechatronics applications: [a pendulum driven by a DC motor](@ref motor-pendulum).
- Electrical and electronic circuit simulation: [a series RLC circuit](@ref rlc-circuit).
- Biomedical systems modeling: [the arterial Windkessel](@ref biomedical-windkessel).
- Chemical process simulation: [coupled tanks](@ref hydraulic-coupled-tanks) and [a chemical reaction network](@ref chemical-reaction-network).

## Installation

To start using the Bond Graph Toolkit, you need to have Julia installed on your system. Follow the official [Julia documentation](https://docs.julialang.org/en/v1/manual/getting-started/) for instructions on installing Julia.

Once Julia is installed, you can add the Bond Graph Toolkit to your project by running the following command in the Julia REPL or the Julia package manager:

```julia
pkg> add https://github.com/felipecastrotc/BondGraphToolkit.jl
```

For more detailed installation instructions, see [Installation](@ref install-guide).

## Getting Started

We recommend checking out the Getting Started to quickly get started with Bond Graph Toolkit. It will guide you through library installation and the modeling and simulation of a simple dynamic system. It provides step-by-step instructions and code examples to help you grasp the basics of using the library.

## Documentation Sections

- [User guide](@ref concepts): the conventions behind bonds, elements and junctions, and how simulation handles causality.
- Examples: classical dynamic systems in the mechanical and electrical domains, multiphysics systems, and one page per feature.
- [Library - Public](@ref public-documentation): Detailed documentation of all public functions provided by the Bond Graph Toolkit.
- [Library - Internals](@ref internal-documentation): Detailed documentation of all internal functions and types used by the Bond Graph Toolkit.

## Acknowledgements

BondGraphToolkit.jl grew out of the author's doctoral research at the [School of Mechanical Engineering (FEM)](https://www.fem.unicamp.br/) of the [University of Campinas (UNICAMP)](https://www.unicamp.br/), supervised by [Prof. Alberto Luiz Serpa](https://bv.fapesp.br/en/pesquisador/4509/), and carried out at [CEPETRO](https://www.cepetro.unicamp.br/), UNICAMP's Center for Energy and Petroleum Studies. The doctorate was supported by the São Paulo Research Foundation (FAPESP), [process 2019/14597-5](https://bv.fapesp.br/en/bolsas/187752/modelling-the-dynamical-behavior-in-electrical-submersible-pump-esp-for-liquid-liquid-flow/).

Opinions, hypotheses and conclusions or recommendations expressed in this material are the responsibility of the author and do not necessarily reflect the views of FAPESP.

## Library

```@contents
Pages = ["lib/public.md", "lib/internals.md"]
```

### [Index](@id main-index)

```@index
Pages = ["lib/public.md"]
```