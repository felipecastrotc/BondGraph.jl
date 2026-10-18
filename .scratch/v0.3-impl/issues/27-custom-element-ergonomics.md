# 27: Simpler custom elements: extra connectors, hand-written elements, submodels

**What to build:** since ADR 0002, a hand-written element (one with an extra connector, like the adapters) needs two power ports and three tie equations, because `bond_ports` is internal. Proposed, in order:

1. Export `bond_ports` (it removes the tie equations from hand-written elements).
2. Let `@oneport`/`@twoport` accept extra connectors (a flange, a pin), so the adapters themselves could be written with `@oneport`.
3. A submodel macro only if ticket 26's example shows it saves something: a junction is already a reusable submodel with a free port (`Junction1(:power => :in, r, c; name)`), and multi-junction submodels are a plain MTK `System`.

Multiport fields are ticket 25.

**Blocked by:** None.

**Status:** needs-info

Registered 2026-10-18 after the usability recheck (`.scratch/api-usability/report.md`, recheck section). Maintainer: move macro improvements to a later version, after 0.3. Open question: which custom elements are needed (one-ports with extra connectors, or reusable submodels)?

## Comments
