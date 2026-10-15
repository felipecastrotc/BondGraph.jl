# Karnopp, Margolis & Rosenberg — System Dynamics: Modeling, Simulation, and Control of Mechatronic Systems, 5th ed., John Wiley & Sons 2012: Ch. 3

PDF↔book page offset: book p. N = PDF p. N+10 (e.g. book p. 40 = PDF p. 50; the PDF's TOC page "v" is PDF p. 4).

Notation: the book writes ports as `e/f` with a half-arrow; `→` below means a bond whose half-arrow points in that direction. Causal stroke is written `|` at one end of a bond. Rule (confirmed by Table 3.6/3.7 and Fig. 3.13): the effort signal travels TOWARD the stroke, so the stroke sits at the end of the bond where effort is an INPUT, and flow travels away from the stroke.

## Topics covered
- §3.1 Basic 1-ports: R (resistor), C (capacitor/compliance), I (inertia), sources Se/Sf; energy of C and I (pp. 37-50, Figs 3.1-3.7, Tables 3.1-3.4)
- §3.2 Basic 2-ports: transformer TF, gyrator GY, modulated MTF/MGY, displacement-modulated TF (pp. 50-57, Figs 3.8-3.10)
- §3.3 3-port junctions: 0-junction, 1-junction, n-port generalization, 2-port junctions as sign flippers (pp. 57-63, Figs 3.11-3.12, Table 3.5)
- §3.4 Causality for 1-, 2-, 3-ports, integral vs derivative (pp. 63-67, Tables 3.6-3.7)
- §3.5 Causality and block diagrams (pp. 67-71, Figs 3.13-3.16)
- Problems (pp. 71-76): only those showing a convention are cited (3-1, 3-6, 3-12, 3-15, 3-18)

## Definitions
- 1-port: element addressed through one power port, one (e, f) pair. (§3.1, p. 37)
- Bond graph elements are idealizations of physical effects, not necessarily one-to-one with components. (Ch. 3 intro, p. 37)
- 1-port resistor R: a static function relates e and f. (§3.1, p. 38, Fig. 3.1)
- 1-port capacitor C (compliance): static relation between effort e and displacement q, with q = ∫f dt. (§3.1, p. 40, Fig. 3.2)
- 1-port inertia I: static relation between momentum p and flow f, with p = ∫e dt. (§3.1, p. 43, Fig. 3.4)
- Variable naming: p = momentum, q = displacement; dq = f dt, dp = e dt (referred to as Eqs 2.2a/2.2b of Ch. 2, not in range). (§3.1, pp. 43-44)
- Tetrahedron of state: R relates e-f, C relates e-q, I relates f-p; the p-q edge would be the memristor, which is not a basic element. (§3.1, p. 48, Fig. 3.6 and footnote)
- Effort source Se: e(t) given, f arbitrary. Flow source Sf: f(t) given, e arbitrary. Domain-named variants exist (S_F, S_V, S_P, S_Q, S_τ, S_ω, S_i); book recommends the generic Se/Sf; SE/SF in computer programs. (§3.1, p. 49, Table 3.4)
- Transformer TF: ideal 2-port, power conserving, same power variable types on both sides (or a transducer when crossing domains). (§3.2, pp. 51-52, Fig. 3.8)
- Gyrator GY: ideal 2-port, power conserving, relates effort at one port to flow at the other. (§3.2, p. 53, Fig. 3.9)
- Modulated TF / GY (MTF, MGY): modulus m or r is a signal on an activated bond; no power is associated with the signal bond. (§3.2, p. 55)
- 0-junction (flow junction, common-effort junction, "parallel"): all efforts equal. (§3.3, p. 58)
- 1-junction (effort junction, common-flow junction, "series"): all flows equal. (§3.3, p. 60)
- Junctions are ideal: power neither dissipated nor stored. (§3.3, p. 58)
- Dual elements: identical type of constitutive law with e and f roles swapped (C/I; 0/1). (§3.4.1, p. 65 footnote)
- Integral causality / derivative causality for C and I: see Causality rules. (§3.4.1, p. 65)

## Sign conventions (power direction, half-arrows, signs at junctions/elements and through TF/GY)
- Half-arrow shows the positive direction of power: power flows in the arrow direction when the product ef is positive. (§3.1, p. 39)
- R: choose the half-arrow pointing TOWARD the resistor ("arbitrary but useful rule" for passive elements). Then R, G > 0 and the law lies in the 1st/3rd quadrant of e-f; it passes through the origin. (§3.1, pp. 39-40, Fig. 3.1b, Table 3.1)
- C and I: same convention, half-arrow toward the element, so ef is power flowing into it and E(t)=∫ e f dt + E0 (Eq. 3.1). (§3.1, pp. 42, 44, Eq. 3.1)
- Sources: half-arrow points away from the source so that e f > 0 means the source supplies power. Book says later it may be convenient to point the arrow the opposite way. (§3.1, p. 49, Table 3.4)
- The book does not state in Ch. 3 what each 1-port law becomes if its half-arrow points away from it. It only notes the sign convention is a choice and shows the effect on junction equations. **Inference:** with an outward arrow on an R/C/I the power into the element is -ef and the laws as written in Tables 3.1-3.3 would need a sign (e.g. e = -R f); the book never writes this. (§3.1, p. 39; flagged)
- Half-arrow direction is a sign-bookkeeping reference: for 0-junction, "half-arrow directions have no effect on the strict equality of the efforts"; for 1-junction, same for flows. They only change the signs in the summation equation. (§3.3, pp. 59, 61, Eq. 3.17, 3.19)
- Junction summation rule: plus sign for bonds with inward-pointing half-arrow, minus sign for outward-pointing. (§3.3, p. 59, Eqs 3.16a, 3.18a)
- 0-junction, all three arrows inward: `e1 f1 + e2 f2 + e3 f3 = 0` (3.16); `e1 = e2 = e3` (3.17); `f1 + f2 + f3 = 0` (3.18). (§3.3, p. 58)
- 0-junction, bond 3 outward: `e1 f1 + e2 f2 - e3 f3 = 0` (3.16a); `f1 + f2 - f3 = 0` (3.18a). (§3.3, p. 59)
- 1-junction, all inward: `f1 = f2 = f3` (3.19); `e1 + e2 + e3 = 0` (3.20). Bond 3 outward: `e1 + e2 - e3 = 0` (3.20a). (§3.3, p. 60)
- Example with bond 1 inward, bond 2 and bond 3 outward at a 1-junction: `f1 = f2 = f3`, `e1 - e2 - e3 = 0` (3.21). (§3.3, p. 61, Eq. 3.21; note the text introducing it says "outward on bonds 2 and 3" in one place and "inward on 2 and 3 but outward on 3" in another; the equation corresponds to 1 inward, 2 and 3 outward)
- Statement of the n-port rule: common effort (0) or flow (1) on all bonds, algebraic sum of the complementary variable is zero. (§3.3, pp. 62-63)
- 2-port junctions: `→0→` and `→1→` (through directions) are equivalent to a plain bond. With opposed arrows they flip a sign: `e1/f1 →0← e2/f2` gives `e1 = e2, f1 = -f2`; `→1←` gives `f1 = f2, e1 = -e2`. Used when two multiports are joined by one bond but defined with incompatible signs (e.g. spring-to-mass: `C ← 1 → I`). (§3.3, p. 63)
- TF and GY use a "through" power convention: both half-arrows point along the 1→2 direction, so power entering port 1 equals power leaving port 2. Eq. 3.8: `e1 f1 = e2 f2`. (§3.2, p. 51, Eq. 3.8, Fig. 3.8a, Fig. 3.9a)
- **Inference:** the book gives TF/GY laws only for the through convention (arrow 1→ element →arrow 2). Reversed arrows on a TF/GY are not treated in Ch. 3, so sign changes in that case are not stated. (§3.2, pp. 51-53; flagged)
- Port numbering matters: in TF/GY the subscript 1 and 2 correspond to the two ports as drawn in the bond graph symbol; which is "1" is fixed by the symbol `→ TF →` with `e1/f1` on the left. (§3.2, p. 51, Fig. 3.8a)
- Footnote: the TF could equally be defined as `m e1 = e2, f1 = m f2`, which inverts the modulus; the book therefore says the modulus convention is ambiguous. Figure 3.10 note: with this convention the modulus there would be 1/(l cos θ). (§3.2, p. 51 footnote; §3.2, p. 57 footnote)
- The summary table for junctions uses all-inward arrows: `e1 f1 + ...`. Table 3.5: 0: `e1=e2=e3`, `f1+f2+f3=0`; 1: `f1=f2=f3`, `e1+e2+e3=0`. (§3.3, p. 63, Table 3.5)
- Block diagrams: a half-arrow becomes a plus or minus sign at a summation point. (§3.5, p. 68)

## Constitutive equations
- R general: `e = \Phi_R(f)`, `f = \Phi_R^{-1}(e)`; linear: `e = R f`, `f = G e = e/R`, R resistance, G conductance. Examples: `F = bV`, `\tau = c\omega`, `P = RQ`, `e = Ri`. (§3.1, p. 40, Table 3.1)
- Turbulent orifice (nonlinear R): `P1 - P2 = (1/(2A^2)) \rho Q |Q|`, equivalently `Q = A \sqrt{(2/\rho)|P1-P2|} sgn(P1-P2)`; A flow area, rho density. (§3.1, p. 39)
- R energy: dissipated; `ef >= 0` with half-arrow toward R. (§3.1, p. 39)
- C general: `q = \Phi_C(e)`, `e = \Phi_C^{-1}(q)`, with `q = \int f dt`; linear: `q = C e`, `e = q/C`. Spring: `C = 1/k`, `F = kx`, x = ∫V dt; torsional spring `C = 1/k_\tau`; tank, accumulator capacitors are C_T, C_a. (§3.1, pp. 40-42, Table 3.2, Fig. 3.2)
- Isentropic air bladder (nonlinear C): `P = P0 V0^\gamma / V^\gamma`, with V = ∫Q dt; gamma ratio of specific heats. (§3.1, p. 42)
- C energy: `E(t) = \int_0^t e f dt + E_0` (3.1); `E(q) = \int_{q0}^{q} e(q) dq + E_0` (3.2); with E0=0 at e=0: `E(q) = \int_{q0}^q e(q) dq` (3.2a); linear: `E = q^2/(2C)`. The shaded area under e-q is the stored energy, recoverable without loss. (§3.1, pp. 42-44, Eqs 3.1-3.2a, Fig. 3.3)
- I general: `p = \Phi_I(f)`, `f = \Phi_I^{-1}(p)`, with `p = \int e dt` (so `\dot p = e`); linear: `p = I f`, `f = p/I`. Examples: mass `p = mV`, `\dot p = F`; rotation `p_\tau = J\omega`; fluid `p_P = I Q`; inductor `\lambda = L i`. (§3.1, pp. 43-47, Table 3.3, Eqs 3.5-3.6)
- I energy: `E(p) = \int_{p0}^p f(p) dp + E_0` (3.3), `E(p) = \int f(p) dp` (3.3a); linear `E = p^2/(2I)`. Stored energy is recoverable. Potential vs kinetic (mechanical), electric vs magnetic (electrical) naming. (§3.1, pp. 45-46, Eqs 3.3-3.3a, Fig. 3.5)
- Relativistic mass as nonlinear I: `p = \Phi_I(V) = mV/(1 - V^2/c^2)^{1/2}` (3.7). Book prefers generalizing `\lambda = L i` to `\lambda = \Phi_I(i)` over generalizing `e = L di/dt`. (§3.1, p. 47, Eq. 3.7)
- Sources: Se: `e(t)` given, `f(t)` arbitrary; Sf: `f(t)` given, `e(t)` arbitrary. (§3.1, p. 49, Table 3.4)
- Transformer (through convention): `e1 = m e2`, `m f1 = f2` (3.9); m is the modulus. Power check: `e1 f1 = e2 f2` (3.8). Lever: `(b/a) V1 = V2`, `F1 = (b/a) F2`; gear: `(r1/r2) \omega1 = \omega2`, `\tau1 = (r1/r2) \tau2`; hydraulic ram (transducer): `F = A P`, `A V = Q` (3.10), A piston area acts as m. Only one modulus exists; the second law follows from power conservation. (§3.2, pp. 51-52, Eqs 3.9-3.10, Fig. 3.8)
- Gyrator: `e1 = r f2`, `r f1 = e2` (3.11); r is the gyrator modulus; r has units of a resistance but the element is power conserving, not dissipative. Voice coil: `e = T V`, `T i = F` (3.12), T the coupling constant, a single value when power units are consistent (SI recommended). (§3.2, pp. 53-54, Eqs 3.11-3.12, Fig. 3.9)
- Cascades: two gyrators `GY1 (r1) → GY2 (r2)` equal a TF with `e1 = (r1/r2) e3`, `(r1/r2) f1 = f3`; two transformers equal a TF with `e1 = m1 m2 e3`, `m1 m2 f1 = f3`. (§3.2, p. 55)
- Gyrator interchanges effort and flow roles: `GY → I` behaves as `C`, `GY → C` as `I`; `I → GY → I` equals `C`. Hence `I-GY-I2` replaced by C gives an oscillator. (§3.2, p. 55)
- Modulated TF/GY: same laws with m or r a signal; `e1 f1 = e2 f2` still holds exactly, no power on the signal bond. (§3.2, p. 55)
- Displacement-modulated TF, rotating arm: `y = l sin\theta` (3.13), `V_y = (l cos\theta) \omega` (3.14), `(l cos\theta) F = \tau` (3.15). Bond graph: `F/V_y → MTF → \tau/\omega` with modulus signal `l cos\theta`. (§3.2, p. 56, Eqs 3.13-3.15, Fig. 3.10)
- 0-junction: `e1=e2=e3`, `f1+f2+f3=0` (all inward). 1-junction: `f1=f2=f3`, `e1+e2+e3=0` (all inward). (§3.3, p. 63, Table 3.5)
- Mechanical 0-junction example: `V3 = -V1 - V2`, with F the common effort. Electrical: 0 = Kirchhoff current law, 1 = voltage law; mechanical 0 = geometric compatibility of velocities, 1 = force balance (Newton for a mass); hydraulic 0 = flow continuity, 1 = pressure-drop loop. (§3.3, pp. 58, 61)

## Causality rules
- Causal stroke: indicates direction in which the effort signal is oriented; effort goes toward the stroke, flow the other way. Sign-convention half-arrow can be on either end without affecting causality. (§3.4.1, p. 64; Table 3.7 text p. 67)
- Se: only `Se—|` (stroke at far end; effort imposed on the system). Sf: only `Sf |—` (stroke at the source end, flow imposed). Equations `e(t)=E(t)`, `f(t)=F(t)`. (§3.4.1, p. 64, Table 3.6)
- R: either causality. Convention: variable on the left of "=" is the output, the argument is the input. `e = \Phi_R(f)` (stroke at the far end, f in) or `f = \Phi_R^{-1}(e)` (stroke at R, e in). If one of Φ_R or its inverse is multivalued, pick the single-valued form. Linear R is indifferent. (§3.4.1, pp. 64-65, Table 3.6, Fig. 3.13)
- C, integral causality (f input, e output): `e = \Phi_C^{-1}(\int f dt)`. Derivative causality (e input, f output): `f = (d/dt) \Phi_C(e)`. (§3.4.1, p. 65, Eq. 3.22, Table 3.6)
- I, integral causality (e input, f output): `f = \Phi_I^{-1}(\int e dt)`. Derivative causality (f input, e output): `e = (d/dt) \Phi_I(f)`. (§3.4.1, p. 65, Eq. 3.23, Table 3.6)
- The integral/derivative distinction holds for linear C and I too. The book defers the consequences to later chapters but says causality predicts model order before writing equations and aids state-equation setup and simulation. (§3.4.1, p. 65; §3.4.3, p. 67)
- TF: only 2 allowed causalities (not 4): once one of e or f is input, the other is fixed. `⊢TF⊢` : `e1 = m e2`, `f2 = m f1`. `⊣TF⊣` : `f1 = f2/m`, `e2 = e1/m`. (§3.4.2, pp. 65-66, Table 3.7)
- GY: one bond's causality fixes the other. `⊢GY⊣`: `e1 = r f2`, `e2 = r f1`. `⊣GY⊢`: `f1 = e2/r`, `f2 = e1/r`. (Reader note: both bonds in a GY take the same "kind" of stroke position pattern as in Table 3.7; see flags.) (§3.4.2, p. 66, Table 3.7)
- 0-junction: exactly one bond has effort as input to the junction (stroke at the junction end of that bond); all other bonds have effort as output (strokes at far ends). Hence n allowed causalities for an n-port. Example (bond 1 determining): `e2 = e1`, `e3 = e1`, `f1 = -(f2 + f3)` for all-inward arrows. (§3.4.3, pp. 66-67, Table 3.7)
- 1-junction: exactly one bond has flow as input to the junction (effort is output on it, stroke at the far end); all others take effort as inputs. n allowed causalities. Example: `f2 = f1`, `f3 = f1`, `e1 = -(e2 + e3)` with all-inward arrows. (§3.4.3, p. 67, Table 3.7)
- Table 3.7 equations use sign patterns for the arrow layout shown there (all inward); they change with outward arrows by the junction sign rule. **Inference:** signs in these causal assignments follow from Eqs. 3.18a/3.20a by moving the signed term. (§3.4.3, p. 66, Table 3.7)

## Assumptions
- 1-ports are described by a static constitutive function; R, C, I are idealized from real components; linear versions are special cases with a parameter. (§3.1, pp. 37-47)
- Resistive R dissipates (power only flows in): law restricted to quadrants 1 and 3 of e-f. Sources are idealized and can supply unlimited power; not realistic over all ranges (battery example, short-circuit current). (§3.1, pp. 39, 49-50, Fig. 3.7)
- Stored energy in C and I is conserved and recoverable; zero reference chosen where e=0 (C) or f=0 (I). (§3.1, pp. 43, 46)
- TF, GY and junctions are ideal: power conserving, no storage or dissipation; actual devices (rigid lever etc.) approximate this only if massless, rigid and frictionless. (§3.2, p. 51; §3.3, p. 58)
- Modulation (m, r varying in time) preserves power conservation. (§3.2, p. 55)
- Use consistent (SI) units for power in all domains so the two gyrator-coupling constants coincide. (§3.2, p. 54)
- Mass as I is exact for Newtonian mechanics (`p = mV`); electrical inductance linear in flux is a modeling decision. (§3.1, p. 47)

## Worked examples
- Resistors, capacitors, inertias in several domains (circuit element, dashpot, orifice; spring, tank, bladder; inductor, mass, pipe, flywheel): bond graphs `R:R`, `R:b`, `C:C`, `C:1/k`, `C:1/kT`, `I:L`, `I:M`, `I:If`, `I:J`, with variables e/f drawn on the bond. (§3.1, pp. 38, 41, 45, Figs 3.1, 3.2, 3.4)
- Battery with a variable resistor: bond graph `Se --e/i--> R`; ideal source predicts current growing unbounded as R→0, real battery saturates at short-circuit current. (§3.1, pp. 49-50, Fig. 3.7)
- Lever, gear pair, electrical transformer, hydraulic ram: each reduces to `TF` with modulus from kinematics, e.g. lever m = b/a, gear m = r1/r2, ram m = A. (§3.2, pp. 51-52, Fig. 3.8)
- Gyroscope (`GY`), voice coil `e = T V`, `T i = F`. (§3.2, pp. 53-54, Fig. 3.9)
- Arm with rotation: displacement-modulated TF `MTF` with modulus `l cos θ`, relating `F, V_y` to `τ, ω`. (§3.2, p. 56, Fig. 3.10)
- Fig. 3.11: 0-junction in electrical (Kirchhoff node), hydraulic (pipe tee), mechanical (two carts, `V3 = -V1 - V2`) domains; 1-junction in series electrical, hydraulic tee with common flow, and mechanical force balance on a massless cart with common velocity. (§3.3, pp. 58-60, Fig. 3.11)
- Fig. 3.12(a): `Se → 0 → I` with `R` on the 0-junction (parallel RL fed by voltage source e; mechanical analog: force source `S_F` on a damper in parallel with a mass, `V2 = -X13' = V1 - V3`). Fig. 3.12(b): `Sf → 1 → C` with `I` on the 1-junction (series LC with current source, mechanical: velocity source, mass and spring). Equations not written out in the text; reader verifies. (§3.3, pp. 61-62, Fig. 3.12)
- `C ← 1 → I` for spring to mass with a sign flip: both bonds carry the same flow, one effort is the negative of the other. (§3.3, p. 63)
- `I - GY - I2` equals `C`; example showing oscillator with C and I. (§3.2, p. 55)
- Drive-train example: bond graph (engine, flywheels, MTF for gravity, TFs, wheel/axle R, I, C) and the equivalent block diagram. No equations. (§3.5, pp. 69-70, Fig. 3.16)
- Problem 3-1 (absquare dashpot `F = A V |V|`: plot, sign convention, causality): indicates that sign convention is a half-arrow toward the R. (Problems, p. 71)
- Problem 3-6: rigid pipe slug as `I` on a 1-junction with pressures P1, P2, P3 and flow Q2; `p_{P2} = ∫ P2 dt = (ρ L/A) Q2`. (Problems, pp. 72-73)
- Problem 3-15: identity `e/i → GY → (F/V) I = e/i → C` to be verified from Eq. 3.12. (Problems, p. 75)
- Problem 3-18: friction as R, only possible causality for it (sign-of-velocity law), singular at v=0. Statement only, no solution. (Problems, p. 76)

## Flags
- Scanned images: half-arrow and stroke positions in Tables 3.6 and 3.7 are small; I inferred which end carries the stroke from the paired equations and from the text on pp. 64 and 67 ("causal stroke on the end of bond 1 nearest the 0 indicates e1 is an input"). Verify visually on pp. 64 and 66 before relying on stroke placement.
- Eq. 3.21 context (p. 61): the prose around it is slightly inconsistent about which bonds are outward; I reported the equation as printed.
- The book does not state in Ch. 3 what an element's law becomes if its half-arrow points away from it (R, C, I), or what TF/GY laws become for reversed arrows. Those statements in the Sign conventions section marked Inference are mine.
- No Ch. 3 statement about general (non-ideal) multiport R/C/I (deferred to Ch. 7) and no explicit discussion of what happens at a 1-port when causality is inconsistent (conflicts, causal loops); deferred to later chapters.
- Eqs. 2.2a, 2.2b, Fig. 2.2 (tetrahedron of state) are in Ch. 2, outside the range; only their use here was read.
- Table 3.7 shows an acausal graph with port 3 oriented upward on the junctions; the arrows there are all inward, as used in the given equations.
- Mechanical-translation C and I use V and F; nothing in the range addresses signs of displacement variables between bonds.
- The table of Table 3.1 SI units and the minor details of Fig. 3.7 were not transcribed beyond what is relevant.
