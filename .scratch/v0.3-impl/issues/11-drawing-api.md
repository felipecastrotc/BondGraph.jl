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

### Review of v0.3 (`/code-review develop`, 2026-10-17): decisions

1. `REPL` without `[compat]`: kept as is. AutoMerge exempts stdlibs (`_AUTOMERGE_REQUIRE_STDLIB_COMPAT = false` in RegistryCI), and `julia = "1.10"` already bounds a bundled stdlib. AGENTS.md and the spec reworded to exempt bundled stdlibs.
2. "MTK internals" in `src/graph.jl`: not a breach. Checked with `Base.isexported`/`Base.ispublic` on MTK 11.45: `Connection` and `getmetadata` are exported; `get_eqs`, `get_systems`, `get_gui_metadata`, `iscomplete`, `t_nounits`, `D_nounits` are declared `public`. The one non-public import was `getname` (through MTK); it now comes from SymbolicIndexingInterface, which exports it (new direct dependency, already installed with MTK). Ticket 08's note that `Connection` is a near-internal was wrong: it is exported.
3. No `@oneport` element tested against a built-in: added "a custom @oneport matches the built-in Compliance" to `test/test_core.jl` (RC with a hand-written linear C vs `Compliance`).
4. Two-inertia case: added "Two inertias driven by an effort share one flow (DAE)" to `test/test_mechanical.jl` (Se, I = 1 and 2, Rodas5P, closed form f = Et/3); stale "use a DAE solver when one is available" comment removed. State count not asserted (MTK's choice).
5. Terminal viewer: the user checks `bgplot(motor)` in a plain terminal REPL by hand.
6. `De`/`Df` docstrings now give the drawing (`[-1, de]`, a package choice: no power, so the direction doesn't matter); `docs/literature/README.md` no longer says De/Df are unimplemented.
7. Adapters stay plain `@component`s (the macros only add signal connectors); recorded in `src/adapters.jl` and ticket 06.
8. Ticket 05's closed-loop criterion reworded to what is tested (`De`, closed form).
9. `public.md` links to Latexify's docs for the re-exported `latexify`.
10. `version = "0.2.1"`: bumped on `release/0.3.0`, no action now.
11. Junction members accept `[1, x]` (explicit form of `x`): kept.
12. `PowerPort` `guess = 0.0` on `e` and `f`: removed. Not needed: with both gone the suite passes with no initialisation warnings (my earlier claim that DAE initialisation needed them was wrong). The adapters' `s`/`phi` guesses stay.
13. `.scratch/codex/*/report.md`: kept as history.
14. Plain `connect` to `p1` and the `PowerPort` wording: moved to ticket 12 for the user to check with examples. Item 5 (terminal viewer) too.
15. `bondgraph` docstring now says a subsystem with a `power` port is drawn as one element (how hand-written elements are recognised).
16. Translational adapters: pure-MTKSL reference comparisons added next to the closed forms (constant force on a 2 kg `Mass`; velocity for `TranslationalAdapter`, position for `TranslationalPositionAdapter`). 13/13 adapter tests; each new comparison fails when the reference force is perturbed.
17. GLOSSARY: "capacitor" is avoided as the element's name, fine for the physical part; the RC example is unchanged.
18. `const Node` for the node `NamedTuple` type in `src/graph.jl`; other duplication left inline.
19. Kept as is.
20. Kept as is.
21. Internal names in `src/graph.jl` (`source`, `graph`) and the naming standard in general: deferred, the user wants a deeper discussion (ticket 12).
22. `src/graph.jl` stays one file.
23. Ternary-with-control-flow in the layout loop and the `walk! = begin … end` closure rewritten as plain `if` / `function`.
24. Shared test helpers and solver `using`s moved to the top of `test/runtests.jl`, so test files don't depend on each other's order.
