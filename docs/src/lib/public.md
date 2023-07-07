# [Public Documentation](@id public-documentation)


Documentation for `BondGraph.jl`'s public interface.

See the [Internals](@ref internal-documentation) section of the manual covering the internal functions.

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
BondGraph
```

### Bond graph elements

```@docs
Se
Sf
Junction1
Junction0
mGY
mTF
Mass
Spring
Damper
```

### Utilities function

```@docs
generate_graph
simplifysys
isindependent
```