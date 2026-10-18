# [Public Documentation](@id public-documentation)

Documentation for `BondGraphToolkit.jl`'s public interface. `using BondGraphToolkit` also re-exports ModelingToolkit, and `latexify` from [Latexify](https://korsbo.github.io/Latexify.jl/stable/).

See the [Internals](@ref internal-documentation) section for the internal functions.

## Contents

```@contents
Pages = ["public.md"]
Depth = 3
```

## Index

```@index
Pages = ["public.md"]
```

## Public Interface

```@docs
BondGraphToolkit
```

### Power port

```@docs
PowerPort
```

### One-port elements

```@docs
Se
Sf
Resistance
Compliance
Inertia
Mass
Spring
Damper
De
Df
```

### Junctions

```@docs
Junction0
Junction1
```

### Two-ports

```@docs
TF
GY
```

### Custom elements

```@docs
@oneport
@twoport
```

### MTK Standard Library adapters

```@docs
ElectricalAdapter
TranslationalAdapter
TranslationalPositionAdapter
RotationalAdapter
HydraulicAdapter
```

### Equations and drawing

```@docs
bgequations
BondGraph
bgplot
bgsave
```
