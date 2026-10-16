# 09: TikZ export

**What to build:** An author exports a model's bond graph as a TikZ `.tex` file in the style of the docs' figures, positioned from Graphviz's layout, ready for a thesis or paper.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 08.

**Status:** ready-for-agent

- [ ] A TikZ writer uses Graphviz's node positions and draws half-arrows, `0`/`1` junctions and `type:parameter` labels in the docs' figure style.
- [ ] Writing the file needs no LaTeX. Viewing it does.
- [ ] Check: the text of a known model has one path per bond with the half-arrow style, plus the expected nodes.
- [ ] A docs example.
