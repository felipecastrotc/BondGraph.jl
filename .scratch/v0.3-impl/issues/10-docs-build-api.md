# 10: Docs build green and API pages complete

**What to build:** A reader browses a docs site that builds cleanly and documents every v0.3 export with a docstring. Getting started uses the v0.3 API, and every element, port, macro and adapter has a runnable example.

Spec: [../../mtk-upgrade/spec.md](../../mtk-upgrade/spec.md). Branch: `feature/mtk-upgrade`.

**Blocked by:** 02, 03, 04, 05, 06, 07, 08, 09.

**Status:** ready-for-agent

- [ ] The public API page lists every export with its docstring; references to removed names are gone.
- [ ] Getting started and the quick example use the v0.3 API.
- [ ] A causality note: MTK chooses the states and reduces the index; a storage element in derivative causality gives a dummy-derivative variable; a causal conflict raises "structurally singular". Show both on small models.
- [ ] Every element, port, macro and adapter has a runnable docs example (checklist against the export list).
- [ ] The docs build passes locally and in CI.
