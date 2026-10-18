# Confirm the target versions and compat bounds

Type: grilling
Status: resolved

## Question

Confirm the targets recommended by the research note (ModelingToolkit 11, Symbolics 7, SymbolicUtils 4, Julia ≥ 1.10 LTS) and how tight `[compat]` should be (e.g. MTK `"11"` vs `"~11.45"`), given BondGraph relies on documented-but-low-level MTK accessors.

## Answer (grilling, 2026-10-15)

- **Targets:** ModelingToolkit 11 and Julia ≥ 1.10. All of MTK 11.45.4, Symbolics 7 and SymbolicUtils 4 require 1.10. MTK 10 would mean a second migration later.
- **MTK bound:** `ModelingToolkit = "11.45"` (caret, ≥ 11.45.0 and < 12). The lower bound is the tested stack (prototype on 11.45.4 / Julia 1.10.12). We use only the public API (ADR 0001), so no `~` pin is needed.
- **Direct dependencies:** only what the code imports. v0.2's Symbolics, SymbolicUtils, MacroTools and Reexport stay only if the new code needs them, each with a bound tested in the same way. The final list is settled in implementation; Plots/GraphRecipes go with the graph-visualization question.
- **CI:** Julia `1.10` and `1` (latest stable). No lower-bound downgrade job and no nightly.
