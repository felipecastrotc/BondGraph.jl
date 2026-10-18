# Borutzky (ed.), Bond Graph Modelling of Engineering Systems: Theory, Applications and Software Support, Springer 2011: Ch. 3 (Djeziri et al.) and Ch. 4 (Borutzky)

PDF↔book page offset: book p. N = PDF p. N+18 (checked: book p. 107 = PDF 125; book p. 135 = PDF 153). Year 2011 from copyright page (DOI 10.1007/978-1-4419-9368-7). Book pages 105-178; book pp. 177-178 are the Part II intro.

## Topics covered
- Ch. 3 §3.2: LFT form and uncertainty attached to elements R, RS, TF, GY, I, C via modulated sources plus virtual detectors (pp. 106-112). (Ch. 3, §3.2, pp. 106-112)
- Ch. 3 §3.2.4: oleopneumatic quarter-vehicle suspension in LFT form, state-space matrices (pp. 113-115). (Ch. 3, §3.2.4, pp. 113-115)
- Ch. 3 §3.3: ARRs, dualised detectors, uncertain ARRs, adaptive thresholds, sensitivity and detectability indices (pp. 116-123). (Ch. 3, §3.3, pp. 116-123)
- Ch. 3 §3.4 and later (pp. 123-134): mechatronic test bench application, FDI results. Skipped: FDI application. (Ch. 3, §3.4, p. 123)
- Ch. 4 §4.2: incremental bond graph (incBG) models of C, R, TF, modulated source, nonlinear R (pp. 137-141). (Ch. 4, §4.2, pp. 137-141)
- Ch. 4 §4.2.2: output sensitivity from incBG state space (pp. 141-142). (Ch. 4, §4.2.2, pp. 141-142)
- Ch. 4 §4.3: direct and inverse models (pp. 142-143). (Ch. 4, §4.3, pp. 142-143)
- Ch. 4 §4.4: transfer-function parameter sensitivities, examples: coupled tanks, DC motor, I-fields and C-fields (pp. 144-156). (Ch. 4, §4.4, pp. 144-156)
- Ch. 4 §4.5: sensitivities for inverse models via bicausality, examples (pp. 157-166). (Ch. 4, §4.5, pp. 157-166)
- Ch. 4 §4.6: sensitivities of ARR residuals, analog integrator example (pp. 166-172). (Ch. 4, §4.6, pp. 166-172)
- Ch. 4 Appendix: state equations with differential-causality stores (pp. 173-174). (Ch. 4, Appendix, pp. 173-174)

## Definitions
- LFT: nominal part M plus a diagonal uncertainty block Δ via auxiliary signals w, z; linear case x' = Ax + B1 w + B2 u, z = C1 x + D11 w + D12 u, y = C2 x + D21 w + D22 u, w = Δ z (Eq. 3.1). (Ch. 3, §3.2.1, p. 106, Eq. 3.1)
- Parameter uncertainty: additive θ = θn + Δθ, or multiplicative θ = θn(1 + δθ), δθ = Δθ/θn; for the inverse parameter, 1/θ = (1/θn)(1 + δ_{1/θ}) with δ_{1/θ} = −Δθ/(θn + Δθ). (Ch. 3, §3.2.2, p. 107, Eqs. 3.2-3.3)
- Detector De / Df on a bond; "dualising" a detector turns it into a source SSe = D̃e / SSf = D̃f imposing the measured value, used to reach derivative causality. (Ch. 3, §3.3, p. 116, Figs. 3.11-3.12)
- ARR: junction balance with unknowns eliminated, leaving known variables (sources, detectors) and parameters; its value is the residual. (Ch. 3, §3.3, p. 116; Ch. 4, §4.6.1, p. 166)
- Incremental bond graph: same structure as the original BG with every element replaced by its first-order incremental model; bonds carry Δe, Δf (variations), not sensitivities. Power variable v(t) = vn(t) + Δv(t) (Eq. 4.1). (Ch. 4, §4.2, p. 137, Eq. 4.1)
- Product Δe·Δf has the dimension of power, hence "true" bond graph, but is only part of the power change. (Ch. 4, §4.2, p. 137)
- Direct model: BG in preferred integral causality, giving state/output from input and known parameters (Eq. 4.16). Inverse model: input u from state, output y and derivatives of y (Eq. 4.17). (Ch. 4, §4.3, pp. 142-143, Eqs. 4.16-4.17)
- Normalized sensitivity index SI_δi = |wi|/a, sum over i equals 1; fault detectability index DI. (Ch. 3, §3.3.2, pp. 120-121, Eqs. 3.30-3.31)
- Parameter fault Y_i (relative deviation of a parameter, modelled like a multiplicative uncertainty) vs structural fault Y_s (new source, e.g. leak as Sf : Ys). (Ch. 3, §3.3.2.2, pp. 120-122, Figs. 3.14-3.15)
- Source-sensor SS: replaces a source plus its detector at collocated input-output pairs, with bicausality, to build the inverse BG. (Ch. 4, §4.5, pp. 157-161)

## Sign conventions
- Junction equations with half-arrow orientation: sum of b_i·e_i (0-junction flows, 1-junction efforts as in the text) with b_i = ±1 following the half-arrow; written ΣSSf + Σ b_i e_i = 0 and ΣSSe + Σ b_i f_i = 0 for dualised detectors. (Ch. 3, §3.3, p. 117, unnumbered eqs.)
- Uncertain ARR: 0-junction Σ b_i f_in + ΣSf + Σ w_i = 0; 1-junction Σ b_i e_in + ΣSe + Σ w_i = 0. (Ch. 3, §3.3.1.1, p. 118, 3rd step)
- Uncertainty-injecting sources carry a negative sign in the LFT figures: w_R = −δ_R z_R etc. (Eq. 3.18). (Ch. 3, §3.2.4, p. 114, Eq. 3.18)
- Mutual inductance sign rule: with both currents entering at the dotted ends, M12 > 0. (Ch. 4, §4.4.3.1, p. 153)
- Incremental TF: modulated sources MSe: −e2n Δm (1-side) and MSf: f1n Δm (0-side) (Fig. 4.3); sign shown in figure, equations below. (Ch. 4, §4.2.1, p. 140, Fig. 4.3)
- Ch. 3, Fig. 3.4 shows the pipe with flow in and P_in − P_out across a 1-junction (effort difference). (Ch. 3, §3.2.2, p. 108, Fig. 3.4)

## Constitutive equations
- Uncertain R, resistance causality: eR = Rn(1+δR) fR = Rn fR + δR Rn fR = e_Rn + e_Runc (Eq. 3.5); modelled as MSe: wR on a 1-junction beside Rn plus virtual De* : zR. (Ch. 3, §3.2.2, p. 107, Eq. 3.5, Fig. 3.2)
- Uncertain R, conductance causality: fR = (1/Rn)(1+δ_{1/R}) eR = f_Rn + f_Runc (Eq. 3.7), via MSf: w_{1/R} on a 0-side junction. (Ch. 3, §3.2.2, p. 108, Eqs. 3.6-3.7, Fig. 3.3)
- Nonlinear R: nonlinearity lives in the characteristic law, not drawn on the BG. Pipe: ṁ_Rz = (1/Rz) √(Pin − Pout), Rz = 8 ρ_l L_p/(π r_p⁴) (Poiseuille), δ_{1/Rz} = −δ_ρl − δ_Lp + 4 δ_rp; w_Rz = −δ_{1/Rz} √(Pin−Pout)/Rz (Eqs. 3.8-3.11). (Ch. 3, §3.2.2, pp. 108-109, Eqs. 3.8-3.11)
- RS element (active resistance, thermal pseudo-bond-graph; port 1 electrical (u, i), port 2 thermal (T, φ̇); generated heat flow f2RS = e1RS f1RS when both port-1 variables are known): when only one is known, f2RS = Φ(RSn, e1RS)(1+δ_{1/RS}) or Φ(RSn, f1RS)(1+δ_RS), drawn as modulated MSf at the output 0-junction (Eqs. 3.12-3.14). Original text writes f2RS = e1RS·f1RS. (Ch. 3, §3.2.2, pp. 109-110, Eqs. 3.12-3.14)
- Uncertain TF, m causality: e1 = mn(1+δm) e2 = mn(e2 − wb), wb = −δm e2; f2 = mn(1+δm) f1 = mn(f1 + wa), wa = δm f1, δm = Δm/mn (Eq. 3.15). 1/m causality: e2 = (1/mn)(1+δ_{1/m}) e1 = (1/mn)(e1 + wa), wa = δ_{1/m} e1; f1 = (1/mn)(f2 − wb), wb = −δ_{1/m} f2; δ_{1/m} = −Δm/(mn+Δm) (Eq. 3.16). (Ch. 3, §3.2.2, pp. 110-111, Eqs. 3.15-3.16)
- GY, I, C uncertainty: same pattern, tabulated in Fig. 3.8 (r and 1/r causality for GY, derivative and integral causality for I, C). (Ch. 3, §3.2.2, pp. 111-112, Fig. 3.8)
- Incremental C: q = C eC gives ΔeC = (1/Cn) Δq − (ΔC/Cn) e_Cn (Eq. 4.3); MSe on a 1-junction with δC = 1/Cn, zC = e_Cn(t), modulated by ΔC. Exact version: qn+Δq = (Cn+ΔC)(en+Δe). (Ch. 4, §4.2.1, pp. 138-139, Eqs. 4.2-4.4, Fig. 4.1)
- Incremental R: eR = R fR gives ΔeR = Rn Δ fR + ΔR fR (Eq. 4.6) or ΔfR = (1/Rn)(ΔeR − (ΔR/Rn) eR) (Eq. 4.7); source MSe: f_Rn(t)·ΔR on a 1-junction. Same model in either causality (Fig. 4.2). (Ch. 4, §4.2.1, p. 139, Eqs. 4.5-4.7)
- Incremental TF: Δe1 = mn Δe2 + e2n Δm, Δf2 = mn Δf1 + f1n Δm (Eq. 4.8). (Ch. 4, §4.2.1, p. 140, Eq. 4.8)
- Effort-modulated source E(t) = k e(t): ΔE = kn Δe + en Δk, extra MSe: en Δk. (Ch. 4, §4.2.1, p. 140, Eq. 4.9, Fig. 4.4)
- Nonlinear R eR = ΦR(fR, Θ): ΔeR = ∂ΦR/∂fR ΔfR + Σj ∂ΦR/∂Θj ΔΘj (Eq. 4.11); modulated resistor MR with slope ∂Φ/∂f plus one MSe per parameter. Orifice fR = cd A sign(eR) √(2|eR|/ρ), eR > 0: MR : k cd/(2√eR), MSf : (k √eR) Δcd, k := A√(2/ρ). (Ch. 4, §4.2.1, pp. 140-141, Eqs. 4.10-4.12, Figs. 4.5-4.6)
- Sources: incremental source has value zero; 0- and 1-junctions unchanged; other elements differ by extra modulated sources/sinks on junctions. (Ch. 4, §4.2, p. 138)
- Incremental state space: Δẋ = An Δx + B*(Θn) w, Δy = Cn Δx + D*(Θn) w, w = W(t,Θn) ΔΘ with W diagonal (Eqs. 4.13-4.14); A and C same as original BG. L{∂y/∂Θ} = [Cn(sI−An)⁻¹B* + D*] L{W} = F* L{W} (Eq. 4.15). (Ch. 4, §4.2.2, p. 142, Eqs. 4.13-4.15)
- Per-parameter source output wj = δj zj ΔΘj (Eq. 4.21); for C, δ = 1/Cn with z = e_Cn (Fig. 4.1) or δ = 1/Cn² with z = q_n; for R, δR = 1, zR = eR, or δR = 1/Rn, zR = eR in the conductance case. Result: ∂F_iκ/∂Θj = F*_ij δj F_j'κ (Eq. 4.27). (Ch. 4, §4.4, pp. 144-145, Eqs. 4.21-4.27)
- Linear I-field: λ = L i with L = [[L1, M12],[M12, L2]]; Δλ = LΔi + (ΔL)i, incBG: derivative-causality I-field with MSe: (ΔL1) i1', (ΔM12) i2', (ΔL2) i2', (ΔM12) i1' on the two 1-junctions; L{Δi} = (I − sAL)⁻¹ sA L{w} (Eqs. 4.49-4.54). C-field (piezoelectric): [x; q] = [[Cm, dε],[dε, Ce]] [F; u], incBG with MSf: (ΔCm) Ḟ, (Δdε) u̇, etc. (Eqs. 4.61-4.62). (Ch. 4, §4.4.3, pp. 153-156, Figs. 4.16-4.18)
- Inverse model: L u = H* L y, H* = H⁻¹ = C*(sI−A*)⁻¹B*(s) + D*(s) (Eqs. 4.18-4.20). incBG of inverse: ΔẋA*, B1*, B2*, B3* (Eq. 4.65), M1 = H*, M2 (∂W/∂Θi) M3 = ∂H*/∂Θi (Eqs. 4.66-4.71). (Ch. 4, §4.3.2 and §4.5.3, pp. 143, 158-159, Eqs. 4.17-4.20, 4.65-4.71)
- Residual sensitivity: ΔL{r} = F* L{W} ΔΘ, ∂L{ri}/∂Θj = F*_ij L{Wjj} = F*_ij δj Σk F_jk L{uk} (Eqs. 4.95-4.97). (Ch. 4, §4.6.3.1, pp. 168-169, Eqs. 4.95-4.97)
- Adaptive threshold: r = −Σ wi, a = Σ|wi|, −a ≤ r ≤ a; sensitivity of a to |δi| equals |e_in| (1-junction) or |f_in| (0-junction). (Ch. 3, §3.3.1.2 and §3.3.2.1, pp. 118-120, Eqs. 3.26-3.29)

## Causality rules
- Property 1.1: a BG is proper iff it has no dynamic element in derivative causality when in preferred integral causality, and vice versa. (Ch. 3, §3.2.3, p. 111)
- Property 1.2 (structural observability): causal path from every integral-causality I/C to a detector, and all I/C admit derivative causality, using dualised detectors if needed. (Ch. 3, §3.2.3, pp. 111-113)
- LFT element models differ per causality (R resistance vs conductance, TF m vs 1/m, GY r vs 1/r, I/C integral vs derivative). (Ch. 3, §3.2.2, p. 112, Fig. 3.8)
- ARRs are generated from derivative causality by dualising detectors; a C stuck in integral causality creates a causality conflict, so the ARR depends on initial conditions. (Ch. 3, §3.3, pp. 116-117, Figs. 3.13, Eq. 3.22)
- Incremental BG model does not depend on the assignment of causality (R example). (Ch. 4, §4.2.1, p. 139)
- Linear I-field in the incBG works in derivative causality. (Ch. 4, §4.4.3.1, p. 154)
- Inverse BG: find a disjoint set of input-output causal paths, replace sources and detectors by SS, assign/propagate bicausality along the paths, then run SCAP on the rest; structural invertibility if a causal path exists (SISO). (Ch. 4, §4.5, pp. 157-158)
- Inverse models usually give derivative causality on some stores, lowering the order vs the direct model; inputs' time derivatives then enter the integral-causality state equations. (Ch. 4, §4.5.2, p. 158; Appendix, p. 173)

## Assumptions
- LFT modelling requires a proper and observable model. (Ch. 3, §3.2.3, p. 111)
- Sources assumed known without uncertainty, except in closed-loop. (Ch. 3, §3.2.2, p. 107)
- Fault detectability: multiple faults give larger contribution than one (Assumption 3.3). (Ch. 3, §3.3.2.2, p. 122)
- Incremental models use first-order Taylor expansion; resulting incBG is linear even when the BG is nonlinear. Exact (full-variation) incBGs also exist. (Ch. 4, §4.2, p. 138)
- Examples assume linear (or linearised around an operating point) models, zero initial values for Laplace analysis, coenergy variables of integral-causality stores as states. (Ch. 4, §4.2.2 and §4.4, pp. 142, 145)
- Appendix assumes a LTI MIMO BG with 1-port stores only, some in differential causality. (Ch. 4, Appendix, p. 173)

## Worked examples
- Oleopneumatic quarter-vehicle suspension (Csph = ν Pe²/(P0 V0), Pe = Ms g/Sp, ν = 1.4), states x = [p_Mns, p_Ms, q_kp = Δl_kp, q_Csph = V_sph]ᵀ (Eq. 3.17); LFT BG with uncertainties on Csph, Ra, Mns, Ms, Sp; auxiliary eqs. 3.18-3.19, matrices A, B1, C1, D11, Δ in Eqs. 3.20-3.21. (Ch. 3, §3.2.4, pp. 113-115, Figs. 3.9-3.10)
- I-R series 1-junction ARR: Se − I1 dSSf/dt − R1 SSf = 0. (Ch. 3, §3.3, p. 117, Eq. 3.22)
- Coupled hydraulic tanks (Fig. 4.8; C1, R1, C2, R2, MSf Qp, Se 0, detectors p1, p2, Qo): ṗ = A p + B Qp with An = [[−1/(C1R1), 1/(C1R1)],[1/(C2R1), −(1/C2)(1/R1+1/R2)]], Bn = [1/C1; 0], y = [p2; Qo] with Cn = [[0,1],[0,1/R2]] (Eq. 4.28); F = (1/Δ)(1/C1)[1; 1/R2](1/(C2R1)) (Eq. 4.30); ∂F1/∂R2 = −(1/(C1C2R1))(1/Δ²) ∂Δ/∂R2 (Eqs. 4.31-4.32), reproduced as F*1 δR2 F2 from the incBG with MSe: Qo ΔR2 (Eqs. 4.33-4.37). (Ch. 4, §4.4.1, pp. 145-148, Figs. 4.7-4.9)
- Fixed-field DC motor (Se E, I:La, R:Ra, GY:kT, I:Jm, R:Rm, MSe Mload): L ω = F21 L E + F22 L M, F21 = kT/(La Jm Δ), Δ = (s+Ra/La)(s+Rm/Jm) + kT²/(Jm La) (Eqs. 4.39-4.42); incBG gives ΔLω = −(1/(Jm Δ))(s+Ra/La)(Lω)ΔRm = F*21 (Lω) ΔRm, hence ∂F21/∂Rm = F*21 F21 (Eqs. 4.43-4.48); Bode plot for Ra=5 Ω, La=50 mH, kT=0.15 Nm/A, Jm=1e-3 kg m², Rm=1e-3 Nms (Table 4.1). (Ch. 4, §4.4.2, pp. 149-152, Figs. 4.10-4.13)
- Mutually coupled coils (Se E, R1, I-field, R2): ΔL = ΔM12[[0,1],[1,0]] gives ∂F/∂M12 = [[F*12F11+F*11F21],[F*22F11+F*21F21]] ΔM12 L E (Eqs. 4.59-4.60). (Ch. 4, §4.4.3.1, pp. 153-156, Figs. 4.14-4.17)
- Inverse of a linear RC network (analog of the tanks; paths I→fC1→x1→e1 order 1 and E→eR2→fR2→f2 order 0, so inverse order is 1): H* entries 4.73 with d = R1C2 s + 1, ∂H*/∂R2 = [[0,0],[0,−1]] (Eq. 4.74), ∂H*/∂C2 = (s/d²)[[1,−R1],[−R1,R1²]] (Eq. 4.75), reproduced by M2 (∂W/∂Θ) M3 (Eqs. 4.76-4.81). (Ch. 4, §4.5.4, pp. 159-163, Figs. 4.19-4.23)
- Inverse DC motor (SS elements, stores in derivative causality, no states): H* = [[La s+Ra, kT],[−kT, Jm s+Rm]], ∂H*/∂Rm = [[0,0],[0,1]] (Eqs. 4.84-4.89); SISO u1 = E, y1 = ω gives u1 = (RaRm/kT + kT) y1 + (LaRm/kT + RaJm/kT) ẏ1 + (LaJm/kT) ÿ1 (Eq. 4.92). (Ch. 4, §4.5.5, pp. 164-166, Figs. 4.24-4.25)
- ARR of coupled tanks: r2 = k1 sign(p1−p2)√|p1−p2| − C2 ṗ2 − Qo (Eq. 4.94); sensitivity of r2 to C2 is −ṗ2. (Ch. 4, §4.6.1, p. 167, Eqs. 4.93-4.94)
- Analog integrator (op-amp macro-model with MSe : A·Vd, R:Ri, R1, Ro, RL, C): residual r via virtual flow detector Df* : Δr; L∂r/∂Ri = F*11 L i_i, L∂r/∂A = F*12 L Vd, with F*11 = (A/Ro − Cs)(1/a), F*12 = 1/Ro, a = 1 + Ri(1/R1 + Cs) (Eqs. 4.98-4.108). (Ch. 4, §4.6.4, pp. 169-172, Figs. 4.26-4.28)
- Appendix state equations with differential-causality stores: ẋi = Ai xi + B_i1 u + B_i2 u̇, xd = Ad xi + Bd u; combined ẋ = A x + B1 u + B2 u̇, y = C x + D u + D1 u̇ (Eqs. 4.109-4.119). (Ch. 4, Appendix, pp. 173-174)

## Coupled hydraulic tanks in detail (used by `examples/hydraulic_coupled_tanks.jl`)
Read: Ch. 4 pp. 141 (Eq. 4.12, Fig. 4.6), 144-149 (§4.4, §4.4.1, Figs. 4.7-4.9, Eqs. 4.21-4.38), 166-167 (§4.6.1, Eqs. 4.93-4.94). Fig. 4.8 was read from a rendered image of PDF p. 164.
- Plant (Fig. 4.7): pump inflow Qp(t) into tank 1; tank levels h1, h2 and bottom pressures p1, p2; valve 1 between the tank bottoms, valve 2 from tank 2 to the outflow Qo(t). (Ch. 4, §4.4.1, p. 145, Fig. 4.7)
- Valve law: "the nonlinear characteristic of the valves is given by Bernoulli's well-known square root law"; the valves are assumed linearised around an operating point so the model is linear and Laplace transform applies. (Ch. 4, §4.4.1, p. 145)
- Orifice law in sign-safe form: fR = cd A sign(eR) √(2|eR|/ρ), cd discharge coefficient, A orifice cross-section, ρ fluid density (Eq. 4.12); in the incremental model k := A√(2/ρ), MR modulus k cd/(2√eR) (the slope ∂fR/∂eR) and MSf : (k√eR) Δcd, assuming eR > 0 (Fig. 4.6). (Ch. 4, §4.2.1, p. 141, Eq. 4.12, Fig. 4.6)
- Bond graph (Fig. 4.8), half-arrows as drawn: MSf → 0 (bond Qp); 0 → C:C1 (p1, QC1); 0 → 1; 1 → R:R1() (QR1); 1 → 0; 0 → C:C2 (p2, QC2); 0 → 1; 1 → R:R2() (QR2); 1 → Se:0 (Qo). Every bond points downstream, and the atmosphere is an effort source of zero whose bond points into it. Detectors De:p1, De:p2 on the two 0-junctions and Df:Qo on the right 1-junction; the text says the two tank pressures and the flow through the second valve are measured. The "()" after R1, R2 is printed as is, not explained on these pages. (Ch. 4, §4.4.1, p. 146, Fig. 4.8)
- States are the co-energy variables of the stores in integral causality (p. 145), here p1, p2; outputs y = [p2; Qo]. Linear model: ṗ = An p + Bn Qp, An = [[−1/(C1R1), 1/(C1R1)], [1/(C2R1), −(1/C2)(1/R1 + 1/R2)]], Bn = [1/C1; 0]; y = Cn p + Dn Qp, Cn = [[0, 1], [0, 1/R2]], Dn = [0]. (Ch. 4, §4.4.1, pp. 145-146, Eq. 4.28)
- Transfer matrix: F = (1/Δ)(1/C1)(1/(C2R1)) [1; 1/R2], Δ = det(sI − An) (Eqs. 4.29-4.30); F1 := L p2/L Qp, ∂F1/∂R2 (Eqs. 4.31-4.32); incremental BG for ΔR2 (Fig. 4.9) with modulated source w = Qo ΔR2, B* = [0; 1/(C2R2)], D* = [0; −1/R2] (Eq. 4.33), same result (Eqs. 4.34-4.37); ∂F1/∂C1 "lengthy" (Eq. 4.38). (Ch. 4, §4.4.1, pp. 146-149)
- ARR on the right 0-junction: 0 = QR1 − QC2 − Qo (Eq. 4.93); with the nonlinear valve 1, r2 = k1 sign(p1 − p2)√|p1 − p2| − C2 ṗ2 − Qo, k1 "a known constant in the constitutive equation of valve 1" (Eq. 4.94). (Ch. 4, §4.6.1, pp. 166-167)
- No numerical parameter values are given for this example on pp. 144-149.
- **Inference:** the tank capacitance in terms of geometry, C = A/(ρg) for an open tank with vertical walls (p = ρgh = ρgV/A), is standard hydrostatics; the pages read use C1, C2 without expressing them in tank dimensions.
- **Inference:** from the slope in Fig. 4.6, the linearised valve resistance at an operating drop Δp* > 0 is Ri = 2√Δp*/(k cd), the R1, R2 of Eq. 4.28 when the linearisation is done this way; the book does not state the operating point or this formula.
- **Inference:** Eq. 4.93 written for the left 0-junction gives Qp − QC1 − QR1 = 0; with QC = C ṗ and the valve law, the two balances are the nonlinear state equations that the example uses.
- Design choice (example, not the book): nonlinear valves (Eq. 4.12) instead of the linearised R1, R2; a constant pump flow Sf instead of MSf : Qp(t); no detectors; parameter values our own. The sign-safe form of Eq. 4.12 and Eq. 4.94 keeps the law defined when p2 > p1 (backflow through valve 1).

## Ideas usable later in the library (sensors/detectors, sensitivity, uncertain parameters)
- Uncertain parameter as modulated element: an element with parameter θ becomes the nominal element plus a modulated source (MSe or MSf) on an adjoining junction, driven by a signal from a virtual detector, fed back through δ. (Ch. 3, §3.2.2, pp. 107-108, Figs. 3.2-3.3)
- Sensitivity via incBG: sensitivity of an output to Θj equals (sensitivity-network gain F*ij) · δj · (nominal transfer to the variable zj that controls the modulated source). (Ch. 4, §4.4, pp. 144-145, Eq. 4.27)
- **Inference:** an MTK incBG could be built by differentiating each element's constitutive equation symbolically; MTK's own `Differential`/structural tooling could then replace the Laplace-domain matrix route. Eq. 4.11 (partials of Φ with respect to state variable and parameters) is the generic recipe.
- **Inference:** an `ARR`/residual could be generated by adding a detector element (De/Df) and replacing it by its dual (SSe/SSf) when building the system, mirroring Figs. 3.11-3.12; with MTK this means turning a measured variable into an input.
- Detectors: De and Df are non-energy-carrying sensors on a bond (signal "activated" bonds), and with dualisation they become imposed sources SSe/SSf. (Ch. 3, §3.3, p. 116)
- Residual sensitivity with respect to a parameter can be added at junctions of the incBG by adding Δ-effort or Δ-flow at those junctions. (Ch. 4, §4.6.3.2, p. 169)
- **Inference:** the incBG structure (same graph, extra signal sinks per parameter) suggests a `sensitivity_system(sys, params)` builder; Chapter 4 itself gives no software API for this beyond CAMP-G/MATLAB use.

## Flags
- Ch. 3 pp. 124-134 (mechatronic test bench, backlash model, FDI results) were not read in detail; FDI theory is skipped by design. (Ch. 3, §3.4)
- Fig. 3.8 and the incremental figures show source orientations only as small arrows; signs of MSe/MSf in figures were read from the printed labels, not verified against the half-arrow direction. (Ch. 3, §3.2.2, p. 112; Ch. 4, §4.2.1, pp. 139-141)
- Eq. 4.21 text on p. 144: for a capacitor, δ_j and z_j are given as 1/Cn^j and e^j_Cn (Fig. 4.1) or (1/Cn^j)² and q^j_n; the second alternative reads ambiguously in the scan. (Ch. 4, §4.4, p. 144)
- Eq. 3.12 f2RS = e1RS·f1RS is printed as is; the thermal pseudo-bond-graph variable assignment was not cross-checked. (Ch. 3, §3.2.2, p. 109)
- Eq. 3.15 for the TF is written both as e1 = mn(1+δm)e2 and with fictive inputs; sign of wb = −δm e2 is as printed. (Ch. 3, §3.2.2, p. 110)
- Ch. 4 uses the Gawthrop-type SS notation and bicausality without defining the bicausal stroke symbol in this chapter. (Ch. 4, §4.5, p. 157)
- Book states Eq. 4.15 and 4.34 as the same F*; Eq. 4.68 uses "ΔW" with entries δiΔΘi, while Eq. 4.14 uses a plain diagonal W. Not an error but a notation shift. (Ch. 4, §4.5.3, p. 159, Eqs. 4.14, 4.68)
