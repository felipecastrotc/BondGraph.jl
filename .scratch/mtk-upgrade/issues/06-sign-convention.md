# Fix the sign convention and the meaning of bond direction

Type: grilling
Status: resolved
Blocked by: 03

## Inputs

- [research/06-sign-convention-analysis.md](../research/06-sign-convention-analysis.md) (Codex, source-and-algebra review): v0.2 signs junction sums only, equalities unsigned, matching 20-sim; the six examples match their code-built graphs; the DC motor and 2-DoF TikZ drawings disagree with their code; candidate conventions A/B/C for v0.3 (A = accepted hybrid).
- Prototype finding (ticket 04): N-DoF spring-damper bond direction mirrors the physics.

## Question

What does the direction of a bond mean in v0.3 (half-arrow convention), how does it set the signs in junction sums, and does the convention apply only to the sums or also to the equalities? Resolve the `TODO.md` doubt (the DC motor matches the literature only with the sign disabled in the equality) against a reference such as 20-sim's convention.

## Answers so far (grilling, 2026-10-15)

Principle (maintainer): the tool follows bond-graph conventions and never assumes what the modeler wants. v0.3 keeps v0.2's sign and direction convention unchanged; the MTK upgrade must implement it, not trade it away.

- Q1 Junction rule: shared variable unsigned; sum signed by each bond's half-arrow (+ into the junction). Matches v0.2 and 20-sim; closes the `TODO.md` doubt (the sign belongs in the sum, not the equality).
- Q2 Element laws are written in the variables of the bond as drawn (v0.2: `Se: e = E`, `Sf: f = F`, `Mass: ṗ = e`). No built-in reorientation: a bond drawn out of an R gives a negative resistor, as drawn. Library elements document their conventional drawing. The "flow always into the element" proposal (convention A) is withdrawn.
- Q3 Standard 1-junction law (common flow, efforts summed with half-arrow signs); no extra "relative signs" rule.
- Q4 `Junction0` honors signs exactly like `Junction1`. Consequence: ADR 0001's "0-junction = plain `connect`" must be amended (0-junction becomes a signed component; plain `connect` is the special case where all bonds point into the connected components).
- Q5 One convention everywhere: TF/GY keep v0.2/textbook laws in bond variables, port 1 in, port 2 out (`TF: e₁ = n e₂, f₂ = n f₁`; `GY: e₂ = g f₁, e₁ = g f₂`).
- Q6 Textbooks are the reference; DC motor and 2-DoF TikZ drawings get fixed (docs step).
- Q7 MTK standard-library adapters follow the same convention as any element: documented conventional direction, signed like any bond.
- Q8 N-DoF example uses the textbook orientation `0 → sd`; its regression check maps the coupling springs' `q` sign against the v0.2 reference data.

## Answer

v0.3 keeps v0.2's convention. The books agree with it ([docs/literature/README.md](../../../docs/literature/README.md)):
- half-arrow = positive power direction;
- junctions: shared variable unsigned, sum signed (+ into, − out) at both junction types (Karnopp Eqs. 3.16–3.21, pp. 58–61; Breedveld Eqs. 1.8–1.11);
- one-port laws in bond variables, conventionally drawn into R/C/I and out of sources (Karnopp pp. 39–49);
- TF/GY through-oriented, with the laws Karnopp Eqs. 3.9, 3.11.

What a reversed one-port arrow means is not stated by either book. The law keeps its form in the bond's variables: a design choice, following the maintainer's principle.

Prototype rerun (branch `prototype/connection-mechanism`, 2026-10-15): with one signed junction component for both 0 and 1 and element laws in bond variables, all 9 golden cases pass with v0.2's own signs (max |Δ| ≤ 2.3e-10). Also:
- the reversed-R RC gives `q = e^{+t}`;
- the textbook N-DoF (`0 → sd`, initial `q = −1`) matches the v0.2 data with the coupling springs' `q` negated;
- mGY power balance is 0.

ADR 0001 is amended.

Follow-ups for implementation:
- the MTK standard-library adapters (ticket 11) were verified under the withdrawn convention A, so re-verify their signs under this one;
- fix the DC motor and 2-DoF TikZ drawings (docs step);
- update `test/regression/helpers.jl` with the N-DoF `q` sign map.
