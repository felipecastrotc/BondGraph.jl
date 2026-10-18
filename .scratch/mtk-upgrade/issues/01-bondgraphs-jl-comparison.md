# How does BondGraphs.jl compare today?

Type: research
Status: resolved

## Question

What is BondGraphs.jl today (scope, API, how it builds on ModelingToolkit/Catalyst, maintenance status, license, registry name), and where does it overlap with or differ from BondGraphToolkit.jl's goal of an MTK-native dialect with easy custom nonlinear elements and modulated two-ports? Also: what do the General registry naming guidelines and AutoMerge similarity checks imply for names like `BondGraph`, `BondGraphToolkit`, `BondGraphsToolkit` next to `BondGraphs`?

## Answer

See [research/01-bondgraphs-jl-comparison.md](../research/01-bondgraphs-jl-comparison.md) (cited).

- BondGraphs.jl (`BondGraphs`, MIT, Forrest et al., Bioinformatics 2023, v0.5.1) models a bond graph as its own Graphs.jl object and converts to MTK at the end; systems-biology focus; custom elements via a Dict library (self-described "unintuitive"); no GY, no modulated TF/GY; effectively MTK 9 only. Complementary to our MTK-native, engineering/multiphysics dialect rather than a replacement.
- `BondGraph` is not registrable: fails all three AutoMerge similarity checks against `BondGraphs` and collides with its exported `BondGraph` type. `BondGraphToolkit` and other longer candidates pass; soft risk: Python's BondGraphTools may draw a review comment.
