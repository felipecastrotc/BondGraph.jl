# Fix the sign convention and the meaning of bond direction

Type: grilling
Status: open
Blocked by: 03

## Inputs

- [research/06-sign-convention-analysis.md](../research/06-sign-convention-analysis.md) (Codex, source-and-algebra review): v0.2 signs junction sums only, equalities unsigned, matching 20-sim; the six examples match their code-built graphs; the DC motor and 2-DoF TikZ drawings disagree with their code; candidate conventions A/B/C for v0.3 (A = accepted hybrid).
- Prototype finding (ticket 04): N-DoF spring-damper bond direction mirrors the physics.

## Question

What does the direction of a bond mean in v0.3 (half-arrow convention), how does it set the signs in junction sums, and does the convention apply only to the sums or also to the equalities? Resolve the `TODO.md` doubt (the DC motor matches the literature only with the sign disabled in the equality) against a reference such as 20-sim's convention.
