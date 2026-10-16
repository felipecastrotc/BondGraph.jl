# 08: Bond-graph extractor and Graphviz drawing

**What to build:** A modeller extracts the bond graph from a model (elements, junctions, two-ports, and bonds with their half-arrows) and renders it with Graphviz as a textbook-style bond graph: half-arrow heads, `0`/`1` junction nodes and `type:parameter` labels, laid out left to right.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 04.

**Status:** ready-for-agent

- [ ] Junction and two-port systems record their member signs and kinds as system metadata, so the extractor uses no MTK internals.
- [ ] An exported extractor returns plain data: nodes with type, label and kind, and bonds with from/to nodes following the half-arrow.
- [ ] A DOT writer (no dependency) produces half-arrow heads, `0`/`1` junctions, `type:parameter` labels and a left-to-right layout.
- [ ] A package extension on a Graphviz package renders the DOT (SVG in docs and notebooks). Weak dependency with a bound.
- [ ] Checks: the extractor on a known model returns the expected nodes, kinds, labels and bond directions; the DOT text has one edge per bond with half-arrow attributes; a rendering smoke test runs in the extension.
- [ ] A docs example rendering one of the examples' bond graphs.
