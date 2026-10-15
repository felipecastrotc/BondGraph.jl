# Decide the public API and its breaking changes

Type: grilling
Status: open

## Already decided (ticket 04)

- Junction-to-junction bonds via named free ports + explicit `connect`; `Mass` has momentum state `p`; `GY`/`TF` components replace `mGY(...; coneqs)` and the tuple `@named` form.

## Inputs

- [research/10-api-usage-inventory.md](../research/10-api-usage-inventory.md) (Codex): every export with signature, call sites in examples/tests/docs, usage patterns, mechanical names used for electrical elements (RC uses `Spring`/`Damper`, DC motor `Mass` for the inductor), dead exports (`renamevars`, `addsubsys`, `latexify`), and which symbols depend on removed MTK internals.

- [research/11-coupling-with-mtk-stdlib.md](../research/11-coupling-with-mtk-stdlib.md): adapters, `RealInput` forms, sensors, MTKSL dependency, and the `Power`/`Mass`/`Spring`/`Damper` name clashes with MTKSL.

## Question

What is the v0.3 public API: element names (`Mass`/`Spring`/`Damper` vs I/C/R per `GLOSSARY.md`), the `@named sys, con = mGY(...)` tuple form, the future of `@oneport` for custom elements, whether `generate_graph` (Plots/GraphRecipes) becomes a package extension, and the exports with no local definition (`renamevars`, `addsubsys`, `latexify`)?
