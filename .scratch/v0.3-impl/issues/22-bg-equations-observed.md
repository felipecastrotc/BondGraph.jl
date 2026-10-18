# 22: `bg_equations` and parameters that appear only in observed equations

**What to build:** `bg_equations` lists only the state equations, so a parameter that appears only in an observed variable can't be shown or renamed. Found by the ticket-17 Windkessel: the characteristic impedance `Zc` appears only in the entrance pressure, so the example can't read it against the hand-derived equations.

**Blocked by:** None.

**Status:** ready-for-agent

Open question: should `bg_equations` (optionally) list the observed equations too, as `observed(sys)` does?

## Comments

- 2026-10-18 (triage, user agreed): 0.3.x (additive, non-breaking): an `observed = true` keyword on `bg_equations` that also lists the observed equations, renameable like the state equations. Then use it in the Windkessel example for `Zc`.
