# Confirm the target versions and compat bounds

Type: grilling
Status: open

## Question

Confirm the targets recommended by the research note (ModelingToolkit 11, Symbolics 7, SymbolicUtils 4, Julia ≥ 1.10 LTS) and how tight `[compat]` should be (e.g. MTK `"11"` vs `"~11.45"`), given BondGraph relies on documented-but-low-level MTK accessors.
