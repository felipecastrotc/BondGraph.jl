# 11: Drawing API: `bgplot`, `savebg`, textbook layout

**What to build:** A modeller looks at a model's bond graph the way they look at a solution. `bondgraph(motor)` gives the bond graph as data and shows it as text in every REPL. `bgplot(motor)` draws it, the way `plot(sol)` draws a solution. `savebg("motor.svg", motor)` writes it to a file, the way `savefig` does. The drawing follows the textbook layout: the power path through junctions and two-ports runs straight left to right, and each element sits beside its node.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md) (drawing section, amended by this ticket). Branch: `feature/mtk-upgrade`. Replaces the drawing API of tickets 08–09.

**Blocked by:** 08, 09.

**Status:** done (feature/mtk-upgrade)

## Decisions (with the user, 2026-10-16)

- Default view is text. `bondgraph(sys)` returns a `BondGraph` (fields `nodes`, `bonds`, as before) whose display is a list of its bonds (`Se:V ⇀ 1:je`) in every REPL, including VS Code and notebooks. Graphics are explicit.
- `bgplot(x; format = :svg)`, where `x` is a model or a `BondGraph`. `:svg` is the drawing: shown in the VS Code plot pane, Jupyter, Pluto and the docs; in a plain terminal REPL it opens in the system's default viewer, as Plots opens a GR window. `:tikz` and `:dot` give the TikZ or DOT source, displayed as plain text.
- `savebg(path, x)`, where `x` is a model, a `BondGraph` or a `bgplot` result. The format comes from the extension: `.svg`, `.png`, `.pdf`, `.tex` (TikZ), `.dot`.
- `bondgraph_dot`, `bondgraph_svg` and `bondgraph_tikz` are removed. The names were confusing (user).
- Rejected: `plot(motor)` and `print(motor)`. Both are type piracy on MTK's `System`, and `print` would replace MTK's display of every model. A Plots recipe (`plot(bondgraph(motor))`) is deferred until someone needs the bond graph inside a Plots layout.
- Graphviz_jll becomes a regular dependency instead of an extension (batteries included; the user hit "Package Graphviz_jll not found"). Reverses spec user story 43. Cost: about 8 MB of binaries plus Cairo/Pango, which Plots usually installs anyway.
- Layout: elements take a free horizontal side of their node first (left, then right), then above, then below; further elements go in the next column. A textbook-style layout rather than a force-directed one (neato/fdp were tried and looked worse).

## Acceptance criteria

- [x] `bondgraph(sys)` displays as a text list of bonds; `g.nodes`, `g.bonds` unchanged.
- [x] `bgplot` gives SVG (image display), TikZ and DOT (text display); unknown formats error.
- [ ] In a terminal REPL, displaying an SVG `bgplot` opens the default viewer. (Not yet checked by hand; no automated test, it opens a window.)
- [x] `savebg` writes each format by extension and errors on an unknown one.
- [x] Graphviz_jll is a dependency with a compat bound; the extension is gone; no `using Graphviz_jll` in docs, examples or tests.
- [x] Textbook layout, checked on the DC motor: `Se` left of its junction on the power path, `R` above it.
- [x] Drawing docs page, API page, AGENTS.md and spec updated.

## Comments

- 2026-10-16: done. `BondGraph` and `BondGraphPlot` are unexported types (documented on the internals page). `print(bgplot(...))` writes the source; `showable("image/svg+xml", p)` is true only for `:svg`, so VS Code, Jupyter, Pluto and Documenter show the drawing, and TikZ/DOT show as text. The terminal viewer is a `display(::REPL.REPLDisplay, ::BondGraphPlot)` method that writes a temporary SVG and runs `open` / `xdg-open` / `start`; it has no automated test because it opens a window. REPL is a new stdlib dependency for that method. 144/144 tests pass; docs build.
