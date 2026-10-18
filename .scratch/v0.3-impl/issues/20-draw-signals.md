# 20: Draw signals and modulated elements in `bgplot`

**What to build:** `bgplot` draws power bonds only. A model with signals (a `Df` feeding a controller block that drives a modulated `Se`) loses the signal path in the drawing, and the modulated source is labelled `Se`, not `MSe`. Textbooks draw activated bonds as full arrows and label modulated elements MSe/MSf/MTF/MGY (Karnopp p. 33, pp. 51–55; Breedveld Fig. 1.1).

**Blocked by:** None.

**Status:** needs-info

Found by the ticket-14 example (DC motor PI control), 2026-10-17. Open questions: draw MTKSL blocks as nodes, or only the signal arrows between bond-graph elements? Label `MSe` from whether the element reads `input`?

## Comments

- 2026-10-17: ticket 15 (motor-driven pendulum) hits the same gap (Df → ∫ → MTF not drawn), and the layout overlaps the `I:mx` label with the `my` bond.
- 2026-10-17: ticket 18 (coupled tanks): a custom element is labelled by its constructor name (`Orifice:valve1`), not by its bond-graph type (R); and only the order of junction members steers the layout (`atm` listed first to sit at the right end).
- 2026-10-18 (triage, user agreed): after v0.3. Two design questions open before it can be specified: draw MTKSL blocks as nodes or only the signal arrows, and how to detect a modulated element (reads `input`?). Custom-element labels (`Orifice` → `R`) and the layout overlap belong here too.
