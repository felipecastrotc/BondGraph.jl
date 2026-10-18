# Bond-graph literature notes

Notes on the two references of the maintainer's thesis, for checking BondGraphToolkit.jl's conventions. These are paraphrases with page references. The PDFs live in `literature/`, which is git-ignored and never committed. Page numbers are the **printed book pages**.

Each note marks the reader's own interpretation with **Inference:**. In this index, "design choice" marks a decision of this project, not something the books state.

## Sources

| Note | Source | Range read | PDF offset |
|---|---|---|---|
| [karnopp-ch1-2.md](karnopp-ch1-2.md) | Karnopp, Margolis & Rosenberg, *System Dynamics: Modeling, Simulation, and Control of Mechatronic Systems*, 5th ed., Wiley 2012, Ch. 1–2 | pp. 4–36 | book p. N = PDF p. N+10 |
| [karnopp-ch3.md](karnopp-ch3.md) | same, Ch. 3 *Basic Bond Graph Elements* | pp. 37–76 | N+10 |
| [borutzky-ch1-breedveld.md](borutzky-ch1-breedveld.md) | Breedveld, *Concept-Oriented Modeling of Dynamic Behavior*, in Borutzky (ed.), *Bond Graph Modelling of Engineering Systems*, Springer 2011, Ch. 1 | pp. 3–52 | N+18 |
| [borutzky-ch2-reduction.md](borutzky-ch2-reduction.md) | Louca et al., *Energy-Based Bond Graph Model Reduction*, same book, Ch. 2 (short, low priority) | pp. 53–104 | N+18 |
| [borutzky-ch3-4.md](borutzky-ch3-4.md) | Djeziri et al., *LFT Bond Graph … Fault Detection* (Ch. 3), and Borutzky, *Incremental Bond Graphs* (Ch. 4) | pp. 105–178 | N+18 |
| [gawthrop-crampin-2014.md](gawthrop-crampin-2014.md) | Gawthrop & Crampin, *Energy-based analysis of biochemical cycles using bond graphs*, Proc. R. Soc. A 470: 20140459, 2014 (arXiv:1406.2447 v2, open access) | §1–§3.1, Eqs. 1–35 (arXiv numbering) | n/a |
| [ctms-dc-motor-speed.md](ctms-dc-motor-speed.md) | Michigan CTMS, *DC Motor Speed: System Modeling* and *PID Controller Design* (read as Internet Archive snapshots, 2022-11-19 and 2023-06-10; text and MATLAB code, not the figures) | both pages | n/a |
| [westerhof-2009-windkessel.md](westerhof-2009-windkessel.md) | Westerhof, Lankhaar & Westerhof, *The arterial Windkessel*, Med Biol Eng Comput 47:131–141, 2009 (Springer HTML full text) | §1–§8, Figs. 3, 5, Tables 1–2 | n/a |

## Topic index

| Topic | Where |
|---|---|
| Power, effort, flow, P = e·f | Karnopp §2.2 pp. 20–25 (Eq. 2.1); Breedveld §1.2 |
| Energy variables p = ∫e dt, q = ∫f dt | Karnopp pp. 20–21 (Eqs. 2.2–2.5) |
| Half-arrow meaning | Karnopp §2.2 p. 25, §3.1 p. 39; Breedveld §1.3.3 pp. 15–16; Louca §2.3.2 p. 73 |
| Default orientation of one-ports (into R/C/I, out of sources) | Karnopp §3.1 pp. 39, 42, 44, 49 (not stated by Breedveld) |
| R, C, I laws (linear and nonlinear), energy, co-energy | Karnopp §3.1 pp. 37–50 (Tables 3.1–3.4); Breedveld Table 1.3 p. 18, pp. 38–43 |
| Se, Sf | Karnopp §3.1 p. 49 (Table 3.4) |
| 0/1-junction equations and signs | Karnopp §3.3 pp. 58–61 (Eqs. 3.16–3.21); Breedveld Eqs. 1.8–1.11 pp. 15–16 |
| TF, GY, MTF, MGY | Karnopp §3.2 pp. 51–55 (Eqs. 3.8–3.12); Breedveld Eqs. 1.2–1.7 p. 15, Eqs. 1.17–1.19 |
| Causal stroke, causality rules, integral vs derivative causality | Karnopp §2.4 pp. 30–32, §3.4 pp. 63–67 (Tables 3.6–3.7); Breedveld §1.2.4, §1.4.2–1.4.3 pp. 7–8, 22–25 |
| Activated bonds (signals) | Karnopp p. 33, pp. 51–55; Djeziri Ch. 3 |
| Worked examples: DC motor | Karnopp Fig. 2.3d pp. 25–26; Breedveld Figs. 1.6–1.8 pp. 19–26 |
| Worked examples: junction examples, RL/LC circuits, mass/spring | Karnopp Figs. 3.11–3.12 pp. 58–62 |
| Power/activity diagnostics | Louca §2.2 pp. 55–58 (Eqs. 2.2, 2.5, 2.6) |
| Detectors De/Df, uncertain parameters | Djeziri Ch. 3 §3.2–3.3 pp. 107–118 |
| Parameter sensitivity (incremental bond graphs) | Borutzky Ch. 4 pp. 135–158 |

## Sign convention: what the books state

1. The half-arrow marks the positive direction of power: when e·f > 0, power flows along the arrow. Both ends of a bond share one e and one f. (Karnopp p. 25, p. 27, p. 39)
2. Junctions: the shared variable is never signed. In the sum, + is used for bonds pointing into the junction and − for bonds pointing out of it.
   - 0-junction: `e₁ = e₂ = e₃` and `Σ ±fᵢ = 0`.
   - 1-junction: `f₁ = f₂ = f₃` and `Σ ±eᵢ = 0`.

   "A 0-junction has only a single effort and a 1-junction has only a single flow, independent of the sign convention." (Karnopp Eqs. 3.16–3.21, pp. 58–61, checked visually; Breedveld Eqs. 1.8–1.11)
3. One-port laws are written in the bond's variables, for the usual orientation. The arrow points into R, C and I, so power flowing into them is positive. Sources point away from the source. Karnopp calls this an "arbitrary but useful rule", and adds that it is sometimes convenient to draw a source's arrow the other way. (Karnopp pp. 39, 42, 49)
4. TF and GY use the through orientation (port 1 in, port 2 out): TF `e₁ = m e₂`, `m f₁ = f₂`; GY `e₁ = r f₂`, `r f₁ = e₂`. Reversing port 2 flips one sign in the law. Modulated versions keep the same laws, with m or r supplied by an activated (signal) bond. (Karnopp Eqs. 3.9, 3.11, pp. 51–55; Breedveld Eqs. 1.2–1.7, p. 15)
5. Causality is independent of the half-arrow. (Karnopp §2.4; Breedveld §1.2.4)

**Not stated by either book:** what a one-port law becomes when its arrow points the other way. **Design choice** (ticket 06, Q2): the law keeps its form in the bond's variables, so you get what you drew. For example, an R drawn with its arrow pointing away acts as a negative resistor.

## Code comparison (v0.2, `develop`)

Callers were found with the codebase knowledge graph: every junction equation comes from `adjmtx2eqs`, which is called by `generate_bg_eqs!` from the MTK override `expand_connections` (`src/mtk.jl:33`). All elements are built in `src/bg.jl`.

| Concept | Books | Code | Status |
|---|---|---|---|
| Bond direction in element ports (v0.3, ADR 0002) | one e and one f per bond (item 1); opposed arrows `→0←` give `e₁ = e₂`, `f₁ = −f₂`, and `→1←` give `f₁ = f₂`, `e₁ = −e₂` (Karnopp §3.3 p. 63) | `bond_ports`, `src/BondGraphToolkit.jl`: `e ~ power.e`, `pout.e ~ power.e`, `f ~ power.f - pout.f`, so `pout` is the `→0←` form | **Agrees**. **Design choice:** flipping the flow rather than the effort; the books allow either, and MTK's `connect` (a 0-junction) can only express `→0←` |
| Junction signs: in the sum, unsigned in the equality | item 2 above | `src/utils/graph.jl:198-268`: −1 for a bond leaving, +1 for one arriving, and `s = apply_equalityeqs ? 1 : s` (l. 249) | **Agrees** |
| Bond direction set by the user | half-arrow | `Junction0`/`Junction1` turn the signs into the direction of `connect` (`src/bg.jl:172-176`, `220-224`) | **Agrees** |
| Se, Sf | `e = E(t)`, `f = F(t)` | `src/bg.jl:63`, `:107` | **Agrees** |
| R | `e = R f` | `Damper`, `src/bg.jl:438` | **Agrees** (linear only; nonlinear R needs a custom element) |
| C | `e = q/C`, `q̇ = f` | `Spring`, `src/bg.jl:412` | **Agrees** |
| I | state p: `ṗ = e`, `f = p/I` (nonlinear: `f = Φ_I⁻¹(p)`) | `Mass`: `D(power.f) ~ power.e / I`, `src/bg.jl:383` (state is f) | **Differs**: works for a linear I, but a nonlinear I can't be written in p. The prototype (ticket 04) already moves to p. |
| TF | `e₁ = m e₂`, `f₂ = m f₁` | `mTF`: `f₁ * r ~ f₂, e₂ * r ~ e₁`, `src/bg.jl:331` | **Agrees** |
| GY | `e₁ = r f₂`, `e₂ = r f₁` | `mGY`: `e₂ ~ g * f₁, e₁ ~ g * f₂`, `src/bg.jl:268` | **Agrees** |
| Modulated TF/GY | modulus from a signal | modulus is any symbolic expression | **Agrees** (the expression plays the role of the activated bond) |
| Causality | stroke rules, integral preferred | none; MTK's structural simplification picks the states and reduces the index | **No implementation** (design: left to MTK) |
| Fields (multiport C/I/R) | Breedveld, Karnopp (later chapters) | only one-ports via `@oneport` | **No implementation** |
| Detectors De/Df | Djeziri Ch. 3 | `De`: `f = 0`, outputs `e` on `output::RealOutput`; `Df`: `e = 0`, outputs `f`, `src/elements.jl` | **Agrees** |
| Power/activity diagnostics | Louca Ch. 2 | none | **No implementation** |
| Sensitivity | Borutzky Ch. 4 | none | **No implementation** (MTK/SciML sensitivity covers it later) |
| DC motor drawing | a source's arrow points away from the source (Karnopp p. 49) | the code bonds the torque source τ into the shaft junction (as the books do), but the TikZ draws `jm → τ` | **Differs** (drawing only): fix the TikZ (ticket 06, Q6) |

## Smallest implementation steps (no code changed yet)

Each step names the check that would catch it breaking.

1. **One signed junction component for both 0- and 1-junctions.** Port variables are the bond's variables. The equality is unsigned, the sum is signed, and the signs come from the user's half-arrows. Check: the 9 golden cases, plus a power balance `Σ sₖ eₖ fₖ = 0` at every junction along a simulated trajectory (Karnopp Eq. 3.16/3.16a).
2. **Keep the v0.2 one-port laws** (Se, Sf, R, C). Check: a one-junction model with R drawn the usual way dissipates (`e·f ≥ 0` on its bond), and the same R drawn reversed delivers power. This pins down the design choice above.
3. **I with state p** (`ṗ = e`, `f = p/I`). Check: the energy `p²/2I` satisfies `dE/dt = e·f` along the trajectory, and the golden 1-DoF case still matches.
4. **TF/GY with the v0.2 laws, port 1 in and port 2 out, constant or modulated.** Check: `e₁f₁ − e₂f₂ = 0` along a trajectory. The prototype already measured 0 for the modulated GY.
5. **No causality code.** Check: one model with dependent states (two inertias on one 1-junction, which puts one of them in derivative causality, Karnopp §3.4 pp. 63–65) compiles to a single state with MTK.
6. **Examples and drawings follow the books:** DC motor per Karnopp Fig. 2.3d, N-DoF with `0 → sd` (ticket 06, Q6/Q8). Check: the regression helpers map the N-DoF `q` sign.

Out of scope for v0.3, but recorded: activity diagnostics, sensitivity via incremental bond graphs, multiport fields.
