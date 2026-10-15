# Which MTK 11 dependencies are AGPL, and what does that mean for us?

Type: research
Status: resolved

## Question

Exactly which packages pulled in by ModelingToolkit v11 are AGPL-3.0 (or other copyleft), since which version, and what that implies for a package that depends on full ModelingToolkit (its own license choice, users, a software paper). What does ModelingToolkitBase (MIT) provide and lack for bond-graph models, concretely for index reduction, tearing, alias elimination, and simulating DAEs from junction structures? Is there an official SciML statement or FAQ on this?

## Answer

See [research/02-licensing-facts.md](../research/02-licensing-facts.md) (cited). Key facts re-verified by the parent session (StateSelection LICENSE = AGPL-3.0; MTK v11.45.4 Project.toml hard-depends on StateSelection and ModelingToolkitTearing).

- Full MTK ≥ v11.0.0 hard-depends on two AGPL-3.0 packages (StateSelection, ModelingToolkitTearing); MTK itself and ModelingToolkitBase stay MIT. No other (A)GPL in the closure; some weak copyleft (MPL-2.0, LGPL-2.1+) either way.
- SciML gives no legal guidance; maintainers' personal views diverge (keep your license but users who build/redistribute follow AGPL vs. a dependent library "needs to be AGPL"). Catalyst.jl precedent: depend on ModelingToolkitBase only to stay MIT; users add `using ModelingToolkit` for full reduction, which swaps in the full `mtkcompile` pipeline.
- ModelingToolkitBase alone has no index reduction, alias elimination, or tearing: index-1 junction models simulate (as mass-matrix DAEs); higher-index cases (rigidly coupled inertias, I in derivative causality) fail or misinitialize. Full MTK reduces them to ODEs.
- JOSS only requires an OSI license on the submitted software; Apache-2.0, MIT and AGPL-3.0 all qualify; no rule on dependency licenses.
- Open: no official FAQ; the real BondGraph examples not yet run on ModelingToolkitBase alone.
