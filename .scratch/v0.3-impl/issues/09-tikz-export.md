# 09: TikZ export

**What to build:** An author exports a model's bond graph as a TikZ `.tex` file in the style of the docs' figures, positioned from Graphviz's layout, ready for a thesis or paper.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 08.

**Status:** done (feature/mtk-upgrade)

- [x] A TikZ writer uses Graphviz's node positions and draws half-arrows, `0`/`1` junctions and `type:parameter` labels in the docs' figure style.
- [x] Writing the file needs no LaTeX. Viewing it does.
- [x] Check: the text of a known model has one path per bond with the half-arrow style, plus the expected nodes.
- [x] A docs example.

## Comments

- 2026-10-15: done. `bondgraph_tikz(sys)` (in the Graphviz extension) lays out with `dot -Tplain`, writes a standalone document with the docs' `bgarrow` half-arrow, and puts names in math mode via `latexify` (so `τ` becomes `$\tau$`, as in the docs' `C: \(C\)` labels). Checked by compiling the DC motor with `pdflatex` and viewing it. Docs page for tickets 08-09: `examples/drawing.jl`.
