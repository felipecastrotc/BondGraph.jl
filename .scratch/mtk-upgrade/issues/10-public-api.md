# Decide the public API and its breaking changes

Type: grilling
Status: open

## Already decided (ticket 04)

- Junction-to-junction bonds via named free ports + explicit `connect`; `Mass` has momentum state `p`; `GY`/`TF` components replace `mGY(...; coneqs)` and the tuple `@named` form.

## Question

What is the v0.3 public API: element names (`Mass`/`Spring`/`Damper` vs I/C/R per `GLOSSARY.md`), the `@named sys, con = mGY(...)` tuple form, the future of `@oneport` for custom elements, whether `generate_graph` (Plots/GraphRecipes) becomes a package extension, and the exports with no local definition (`renamevars`, `addsubsys`, `latexify`)?
