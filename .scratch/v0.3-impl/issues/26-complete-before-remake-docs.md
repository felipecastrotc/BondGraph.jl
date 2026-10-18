# 26: Document which symbols work for `remake` (uncompiled `model.x` doesn't)

**What to build:** the examples index solutions through the uncompiled model (`sol[model.mass.f]`), which works. The same style fails for parameter overrides: `remake(prob; p = [parent.A.j.r.R => 4.0])` warns that `parent₊A₊j₊r₊R` is not a parameter and asserts a missing `Initial` parameter. An uncompleted MTK system prefixes its own name; `complete(parent).A.j.r.R` and `sys.A.j.r.R` (compiled) both work. This is upstream MTK behaviour, not a BondGraphToolkit bug. Add the reusable-subsystem example the usability report asks for (U01): a composite with a boundary port, two instances, reading outputs and overriding one instance's parameter with `remake`, using `sys.x` (or `complete`) and saying why.

**Blocked by:** None.

**Status:** needs-triage

Source: `.scratch/api-usability/` finding P05-F02 / U01 (trial 05, `runs/05/attempt_01.jl:37`, `attempt_01.log`); rechecked on the oriented-ports API 2026-10-18 with the same result. P01-F01 (`prob.ps[sym]` for reading a parameter value) fits in the same example.

## Acceptance criteria

- [ ] An example (or a section of an existing one) builds a reusable subsystem with a port and two instances, and remakes one instance's parameter.
- [ ] It states that overrides use the compiled system's (or a `complete`d model's) symbols, and shows `prob.ps[sym]`.
- [ ] Its result is checked against a closed form.

## Comments
