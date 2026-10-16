# [Concepts and Conventions](@id concepts)

BondGraphToolkit follows the conventions of the bond-graph literature, mainly Karnopp, Margolis & Rosenberg, *System Dynamics* (5th ed., Wiley 2012), and Borutzky (ed.), *Bond Graph Modelling of Engineering Systems* (Springer 2011). This page states them once; every element and junction uses them.

## Bonds and power

A bond carries one effort `e` and one flow `f`, and its half-arrow marks the positive direction of power: when `e·f > 0`, power flows along the arrow (Karnopp §2.2). Every element exchanges power through a [`PowerPort`](@ref) with the variables `e` and `f`; its bond to the rest of the model is an ordinary ModelingToolkit `connect`.

## Elements

An element's law is written in the variables of its bond as drawn (Karnopp §3.1):

| Element | Law | Conventional drawing |
|---|---|---|
| `Se(E)` | `e = E` | bond points away from the source |
| `Sf(F)` | `f = F` | bond points away from the source |
| `Resistance(R)` | `e = R f` | bond points into the element |
| `Compliance(C, q)` | `q̇ = f`, `e = q / C` | bond points into the element |
| `Inertia(I, p)` | `ṗ = e`, `f = p / I` | bond points into the element |

The tool never reinterprets a drawing: a resistance whose bond points away from it keeps `e = R f` in that bond's variables, so it delivers power, a negative resistor. The books don't say what a reversed one-port means; this is the package's choice, so that you always get what you drew.

`Mass(m, u)`, `Spring(k, x)` and `Damper(c)` are mechanical shorthands for `Inertia`, `Compliance` and `Resistance`.

## Junctions

At a 0-junction all bonds share one effort and the flows sum to zero; at a 1-junction all bonds share one flow and the efforts sum to zero. In the sum, a bond pointing into the junction counts `+` and one pointing out counts `-`. The shared variable itself is never signed (Karnopp Eqs. 3.16–3.21).

A junction member is written with the sign of its bond:

```julia
@named m = Junction1(se, [-1, R], [-1, C])   # Se ⇀ 1, 1 ⇀ R, 1 ⇀ C
```

A bare member's bond points into the junction (`+1`), and `[-1, x]` means the bond points out of the junction, towards `x`. A `Symbol` creates a free port: `[-1, :out]` on one junction and `:in` on another, joined with `connect(j1.out, j2.in)`, is a bond from the first junction to the second. A `connect` of more than two free ports is itself a 0-junction.

## Two-ports

[`TF`](@ref) and [`GY`](@ref) have port `p1` (bond in) and `p2` (bond out): `TF`: `e₁ = m e₂`, `f₂ = m f₁`; `GY`: `e₁ = r f₂`, `e₂ = r f₁` (Karnopp Eqs. 3.9, 3.11). A modulus that is an expression, or that comes from a signal, makes them modulated (MTF, MGY).

A plain `connect` to an element's `power` port (from a two-port's port or a junction's free port) is a bond pointing into the element, because ModelingToolkit counts every port's flow positive into its component. A source therefore reaches a two-port through a junction, `Se ⇀ 1 ⇀ TF`, as in the books.

## Simulation and causality

A bond-graph model is an ordinary ModelingToolkit system: `mtkcompile` reduces it, and `ODEProblem` and `solve` simulate it. BondGraphToolkit doesn't assign causal strokes. ModelingToolkit chooses the states and reduces the index itself, which covers what causality analysis would tell you about the equations.

**Derivative causality.** Two inertias on one 1-junction share one flow, so only one of them can keep an independent momentum: in bond-graph terms, the other is in derivative causality (Karnopp §3.4). MTK keeps one momentum as a state and adds a "dummy derivative" algebraic variable for the other:

```@example causality
using BondGraphToolkit
@named se = Se(1.0)
@named m1 = Inertia(I = 1.0)
@named m2 = Inertia(I = 2.0)
@named model = Junction1(se, [-1, m1], [-1, m2])
sys = mtkcompile(model)
unknowns(sys)
```

The result is a small DAE (a mass-matrix ODE), so solve it with a DAE-capable solver such as `Rodas5P` from OrdinaryDiffEqRosenbrock rather than `Tsit5`.

**Causal conflicts.** Two flow sources on one 1-junction both try to set the shared flow, a causal conflict. MTK reports it as a structurally singular system:

```@example causality
@named a = Sf(1.0)
@named b = Sf(2.0)
@named conflict = Junction1(a, [-1, b])
try
    mtkcompile(conflict)
catch err
    first(sprint(showerror, err), 120)
end
```

## Signals

Activated bonds carry information but no power. Sources and two-ports without a value read it from `input` (a `RealInput`), and the detectors [`De`](@ref) and [`Df`](@ref) put an effort or a flow on `output` (a `RealOutput`), so MTK Standard Library blocks can drive and read a bond graph. See [Signals and Detectors](@ref signals).
