# Gawthrop & Crampin 2014: Energy-based analysis of biochemical cycles using bond graphs

P. J. Gawthrop and E. J. Crampin, *Energy-based analysis of biochemical cycles using bond graphs*, Proc. R. Soc. A 470: 20140459, 2014. Open access, doi:10.1098/rspa.2014.0459.

Read: the arXiv version, [arXiv:1406.2447v2](https://arxiv.org/abs/1406.2447) (11 Sep 2014), §1–§3.1 (abstract, introduction, Fig. 1, §2.1–2.5, Fig. 2, §3.1; Eqs. 1–35). §3.2 onwards (null spaces, model reduction, the enzyme and switch cycles of §5, software) not read. Section, figure and equation numbers below and page numbers (PDF pages) are the arXiv v2 ones; **not checked** against the journal's typeset version, whose numbering may differ.

Used by: `examples/chemical_reaction_network.jl` (Fig. 1c, A ⇌ B ⇌ C).

## Power variables and bonds (§2.1)

- A bond carries the chemical potential µ (J mol⁻¹) and the molar flow rate v (mol s⁻¹); their product is power, P = µ v (W). Bonds transmit power and do not create, store or dissipate it. The half-arrow gives the direction in which power is counted positive: a sign convention. (§2.1, p. 4)
- Footnote 2: in standard bond-graph terms, chemical potential is the effort (analogous to voltage, force) and molar flow the flow (analogous to current, velocity). (p. 4)

## Species: C component (§2.1)

- Each pool of a species is a C component holding x moles; for A ⇌ B (Fig. 1a), ẋa = −v and ẋb = v (Eq. 1).
- Dilute system in a volume V: µa = µ0a + RT ln(xa/V) (Eq. 2, citing Keener & Sneyd §1.2), rewritten as µa = RT ln(Ka xa) with Ka = e^(µ0a/RT)/V (Eq. 3).
- The C component stores but does not create or dissipate energy. Footnote 3: C stands for "capacitor"; the potential is like a capacitor's voltage, which charges with a net influx. (p. 4, p. 6)
- The paper labels the species component `C` (e.g. `C:A`), not `Ce`.

## Reaction: Re component (§2.1)

- One Re component per reaction relates the reaction flow v to the forward and reverse affinities, Af and Ar (weighted sums of chemical potentials). For A ⇌ B, Af = µa and Ar = µb. (p. 6)
- Marcelin–de Donder formula: v = v⁺ − v⁻, v⁺ = κ e^(Af/RT), v⁻ = κ e^(Ar/RT) (Eq. 4), κ a constant that sets the reaction rate. The de Donder form v⁺/v⁻ = e^(A/RT) with A = Af − Ar (Eq. 5), and the Marcelin form v = κ(e^(Af/RT) − e^(Ar/RT)) (Eq. 6), "used in the sequel".
- The Re component dissipates, but does not create or store, energy. (p. 6)
- Substituting Eq. 3 into Eq. 4 gives first-order mass action: v⁺ = κKa xa, v⁻ = κKb xb, v = κ(Ka xa − Kb xb) (Eq. 7); in rate-constant form v = k⁺xa − k⁻xb (Eq. 8) with k⁺ = κKa, k⁻ = κKb (Eq. 9). The K's are thermodynamic, κ is kinetic; the rate-constant form no longer separates them. (p. 6)
- The paper states (p. 6, after Eq. 7) that "the equilibrium constant is given by Kb/Ka". **Inference:** setting v = 0 in Eq. 7 gives xb/xa = Ka/Kb = k⁺/k⁻, the reciprocal. Either the paper defines the constant the other way round (xa/xb) or this is a slip; the arXiv text doesn't say which. The example uses K₁ = xb/xa = Ka/Kb, derived from Eq. 7, and doesn't quote the paper's sentence.

## Junctions (§2.2)

- 0-junction (common potential): all impinging bonds share one chemical potential; the molar flows sum to zero with + for inward bonds and − for outward ones. Central junction of Fig. 1c: one bond in (v1), two out (vb, v2), so v1 − vb − v2 = 0, i.e. vb = v1 − v2 (Eq. 10). It conserves power (Eq. 11). The left junction gives va = −v1, the right one vc = v2. (p. 7)
- 1-junction (common flow), dual of the 0-junction: flows equal, efforts sum to zero with the same in/out signs. For A + B ⇌ C + D (Fig. 1d), Af = µa + µb (Eq. 12) and Ar = µc + µd (Eq. 13), giving second-order mass action (Eqs. 14–15). (pp. 7–8)
- These signs are the same as Karnopp's (see `karnopp-ch3.md`, §3.3).
- **Inference:** a unimolecular chain such as Fig. 1c needs no 1-junctions; they appear only when a reaction has several reactants or products.

## Bond directions in Fig. 1c (A ⇌ B ⇌ C)

From Eq. 10 and the text around it: each species' bond points out of its 0-junction into the C component (so the species' flow is ẋ, Eq. 30); each reaction's bonds follow the forward reaction, from the reactant's junction into Re and out of Re into the product's junction. Fig. 1c itself was read only as extracted text, so the arrow heads were not checked visually.

## Stoichiometry and other items read

- TF:n represents a stoichiometric coefficient: for A ⇌ 2B (Fig. 1e), Ar = 2µb, vb = 2v (Eq. 17), and v = κ(Ka xa − (Kb xb)²) (Eq. 18). (§2.3)
- An external flow v0 through A ⇌ B ⇌ C (Fig. 1g) gives a non-equilibrium steady state with v1 = v2 = v0. (§2.4)
- Around a cycle (Fig. 1h) the detailed-balance constraint k⁺₁k⁺₂k⁺₃k⁺₄ / (k⁻₁k⁻₂k⁻₃k⁻₄) = 1 (Eq. 19) holds automatically, because each K appears once in the numerator and once in the denominator (Eqs. 20–23). (§2.5)
- Causal strokes give the stoichiometric matrix, Ẋ = N V (Eq. 26); for Fig. 2b (= Fig. 1c), N = [−1 0; 1 −1; 0 1] (Eqs. 28–30). With C components imposing potential on Re and Re imposing flow on the C's, this is integral causality. (§3.1)

## Relevance to BondGraphToolkit

- **Inference:** the paper's conventions are Karnopp's: µ as effort, v as flow, half-arrow as positive power, signed sums at junctions. So the library's `Junction0`/`Junction1` express the chemical bond graphs as drawn, without any change of convention.
- **Inference:** the C component is a nonlinear compliance (displacement x = ∫v dt, effort RT ln(Kx)), expressible with `@oneport`. Re is a nonlinear two-port with equal flows on both bonds and a flow set by both efforts; it isn't a TF or GY, but `@twoport` expresses it (`f1 = κ(e^(e1/RT) − e^(e2/RT))`, `f2 = f1`).
