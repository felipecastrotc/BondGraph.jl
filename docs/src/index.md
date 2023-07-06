# Welcome to Bond Graph Toolkit Documentation

*A bond graph toolkit for Julia.*

The `BondGraphToolkit` is a user-friendly Julia library for modeling dynamic systems with bond graphs. 

**Introduction**: Start with an introduction section that provides an overview of the library and its purpose. Explain what the library does, its main features, and its intended use cases. This section should provide a high-level understanding of the library's capabilities.

## Package Features

The key features of the package include:

1. **Intuitive and Expressive Modeling Syntax**: The package provides an intuitive and expressive syntax that simplifies modeling dynamic systems with bond graphs.
2. **Wide Range of Bond Graph Components**: It supports bond graph components and elements, including inertance, compliance, resistance, gyrators, transformers, and more.
3. **Simulation Capabilities with `DifferentialEquations.jl` and `ModelingToolkit.jl`**: The package seamlessly integrates with `DifferentialEquations.jl` and `ModelingToolkit.jl` packages. Every bond graph component is an `ODESystem`.
4. **Visualize Bond Graph Systems**: It offers a visualization tool that allows to generate visual representations of the bond graph created.
5. **Customization and Modification**: Users can define their specialized bond graph components with a user-friendly API and extend the package's capabilities.
6. **Extensible Architecture for Integration**. Integrating with other libraries and frameworks like `ModelingToolkit.jl` is simple.

## Use Cases

- Control systems design and analysis.
- Robotics and mechatronics applications.
- Electrical and electronic circuit simulation.
- Biomedical systems modeling.
- Chemical process simulation.

## Installation

To start using MyLibrary, you need to have Julia installed on your system. Follow the official [Julia documentation](https://docs.julialang.org/en/v1/manual/getting-started/) for instructions on installing Julia.

Once Julia is installed, you can add the Bond Graph Toolkit to your project by running the following command in the Julia REPL or the Julia package manager:

```julia
pkg> add BondGraphToolkit.jl
```

For more detailed installation instructions, including information about dependencies, please refer to the [Installation](@ref install-guide).

## Getting Started

We recommend checking out the [Getting Started](@ref) to quickly get started with Bond Graph Toolkit. It will guide you through library installation and the modeling and simulation of a simple dynamic system. It provides step-by-step instructions and code examples to help you grasp the basics of using the library.

## Documentation Sections

TODO: Fix the links
- [User guide](@ref): Learn how to model different types of dynamic systems and perform various analyses using MyLibrary.
- [Examples](@ref): Collection of classical examples of dynamic systems in the mechanical and electrical domains and multiphysics systems.
- [Library - Public](@ref public-documentation): Detailed documentation of all public functions provided by the Bond Graph Toolkit.
- [Library - Internals](@ref internal-documentation): Detailed documentation of all internal functions and types used by the Bond Graph Toolkit.

## Versioning and Changelog

For information about the different versions of MyLibrary and the changes introduced in each version, please refer to the [Changelog](@ref).

## Library

```@contents
Pages = ["lib/public.md", "lib/internals.md"]
```

### [Index](@id main-index)

```@index
Pages = ["lib/public.md"]
```