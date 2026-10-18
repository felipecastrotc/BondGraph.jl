# 03: Examples and docs on oriented ports

**What to build:** Convert every example and the user guide to the new syntax and prose, read results from `x.e`, `x.f`. `two_ports.jl` shows `Se ⇀ TF` with a plain `connect`. Close ticket 12's `p1` item.

**Blocked by:** 01, 02.

**Status:** done (feature/mtk-upgrade, uncommitted), user review pending

## Acceptance criteria

- [x] Every example runs with `julia --project=docs`, numbers in the text unchanged (or the change explained).
- [x] The docs build; no `[-1,` left outside the ADRs and the error message.
- [x] Ticket 12's `p1` item checked off with a pointer to ADR 0002.

## Comments

- 2026-10-18 (agent): all examples and the user guide converted (code by the rewriter, prose by hand: each `[-1, x]` / `:a` explanation now names the new form). `two_ports.jl`'s lever is now `connect(F.pout, tf.p1)`, `connect(tf.p2, mass.power)`, no junction; its by-hand text updated. `elements_oneport.jl`'s hand-written element shows the two ports and the `e`, `f` tie equations. All 18 examples run with `julia --project=docs` (none skipped except `tmp_*`); the physics is unchanged, so the numbers in the text are too (the regression data matches). Docs build: exit 0, only the existing large-plot size warning; 114/114 images resolve. `[-1,` remains only in the ADRs, the removed-form error message, its test, and the `tmp_dev_*` scratch.
