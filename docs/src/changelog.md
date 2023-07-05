# [Duffing Equation](@id msd-duffing)

The Duffing equation is an extension of the simple mass-spring-damper system that incorporates a cubic spring. Despite its simple form, the Duffing equation exhibits chaotic behavior, and it is used to model more realistic damped oscillators. The equation is given by:

$m\,\frac{d^2 x}{dt^2} + c\,\frac{d x}{dt} + k\,x + k_3\,x^3 = F(t),$

where $x$ represents the position, $m$ is the mass, $c$ is the damping coefficient, $k$ is the spring stiffness, $k_3$ is the cubic spring stiffness, and $F(t)$ denotes an external force applied to the system.

In standard bond graphs, the capacitance element is linear, while for the Duffing equation we need a non-linear capacitance, cubic spring. Therefore, to model the Duffing equation using the Bond Graph Toolkit, we need to introduce a non-linear capacitance element. The Bond Graph Toolkit allows users to define custom elements that utilize the [`Power`](@ref) connector.

## Bond Graph Toolkit

To begin, we import the Bond Graph Toolkit module for modeling the system and the `DifferentialEquations.jl` package for solving the resulting ordinary differential equation (ODE). Additionally, we import the independent variable `t` from the Bond Graph Toolkit, which will be used to define a custom forcing term.
