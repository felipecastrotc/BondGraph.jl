# Can junctions be expressed with MTK's native connections?

Type: research
Status: open

## Question

On ModelingToolkit v11, can a 0-junction be expressed with MTK's native acausal `connect` semantics (shared effort as the potential, flows summing to zero as a `Flow` variable), and a 1-junction as its dual (shared flow, efforts summing to zero)? How would bond direction / the sign of each bond be expressed (orientation conventions, sign flips in the port, or junction components with signed ports)? What do existing MTK-based bond-graph or port-Hamiltonian efforts (including BondGraphs.jl and ModelingToolkitStandardLibrary) do? List each candidate design with what it costs in ergonomics and what it breaks.
