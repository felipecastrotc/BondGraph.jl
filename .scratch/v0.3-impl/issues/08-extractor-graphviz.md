# 08: Bond-graph extractor and Graphviz drawing

**What to build:** A modeller extracts the bond graph from a model (elements, junctions, two-ports, and bonds with their half-arrows) and renders it with Graphviz as a textbook-style bond graph: half-arrow heads, `0`/`1` junction nodes and `type:parameter` labels, laid out left to right.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 04.

**Status:** done (feature/mtk-upgrade)

- [x] Junction and two-port systems record their member signs and kinds as system metadata, so the extractor uses no MTK internals.
- [x] An exported extractor returns plain data: nodes with type, label and kind, and bonds with from/to nodes following the half-arrow.
- [x] A DOT writer (no dependency) produces half-arrow heads, `0`/`1` junctions, `type:parameter` labels and a left-to-right layout.
- [x] A package extension on a Graphviz package renders the DOT (SVG in docs and notebooks). Weak dependency with a bound.
- [x] Checks: the extractor on a known model returns the expected nodes, kinds, labels and bond directions; the DOT text has one edge per bond with half-arrow attributes; a rendering smoke test runs in the extension.
- [x] A docs example rendering one of the examples' bond graphs.

## Comments

- 2026-10-15: done. `bondgraph(sys)` returns `(; nodes, bonds)`; junctions record `(kind, members)` as system metadata, elements/two-ports are recognised by their `power` / `p1`,`p2` ports and typed by MTK's `@component` GUI metadata. Reading `connect` statements uses MTK's `Connection` type (not an override; the only near-internal). (Correction 2026-10-17: `Connection` is exported by MTK 11, and the accessors used are declared `public`; see ticket 11's review notes.) Labels are `type:name` (the instance name; every R/C/I parameter is called R/C/I, so the name is more informative). `bondgraph_dot` keeps junction/two-port bonds straight (weight) and uses straight lines. Rendering: package extension on `Graphviz_jll` (`bondgraph_svg`). Signals (activated bonds) are not drawn yet. The docs example lands with ticket 09 on one page.
