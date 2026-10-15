# Breedveld, "Concept-Oriented Modeling of Dynamic Behavior", in Borutzky (ed.), Bond Graph Modelling of Engineering Systems: Theory, Applications and Software Support, Springer 2011: Ch. 1
(ISBN 978-1-4419-9367-0, DOI 10.1007/978-1-4419-9368-7; chapter pp. 3-52.)

PDF↔book page offset: book p. N = PDF p. N+18 (PDF p. 20 shows printed "2"; chapter 1 opens on book p. 3 = PDF p. 21). All citations below are BOOK pages.

## Topics covered
- Bond graph as labeled directed graph; bonds, signals, multiports, modulation (§1.2.2, p. 5)
- Ports, effort/flow, power bonds vs pseudo-bonds, signals (§1.2.3, pp. 6-7)
- Causal stroke as bilateral signal flow (§1.2.4, pp. 7-8, Fig. 1.4)
- Domain tables: thermodynamic vs mechanical framework of variables (Tables 1.1, 1.2, pp. 9-10)
- General multiport constitutive form, Eq. (1.1) (p. 10)
- Node categories: storage, source, sensor (§1.3.2, pp. 12-13)
- Resistive, TF, GY, junctions, RS (§1.3.3, pp. 14-17, Eqs. 1.2-1.12)
- Basic elements and Table 1.3 of linear constitutive relations (§1.3.4, pp. 17-18)
- Transforming an iconic model into a bond graph, 8 steps (§1.3.5, p. 19)
- Pendulum + DC motor example, Figs. 1.6-1.8 (pp. 19-26) and Eqs. 1.13-1.19 (pp. 20-21)
- Causal port properties and causality assignment algorithm (§1.4.2-1.4.3, pp. 23-25)
- Modeling vs model transformation, checklist (§1.5.3, §1.5.4, pp. 27-30)
- State, extensive/intensive, equilibrium-determining/establishing variables (§1.5.6, pp. 30-31)
- Port concept, bond graph construction, pendulum reduction with gyristors (§1.5.8-1.5.9, pp. 32-35, Figs. 1.10-1.11)
- Energy, co-energy, Legendre transform, power, multiport C (§1.5.10-1.5.12, pp. 36-43, Eqs. 1.27-1.58)
- Thermodynamic vs mechanical framework, Hamiltonian (§1.5.13, pp. 43-45, Eqs. 1.59-1.63)
- Energy states vs configuration states, symplectic gyrator (§1.5.14, pp. 45-46)
- Conservation/continuity within a domain (§1.5.15, p. 50)

## Definitions
- Bond graph: labeled di-graph with two edge types, bonds and signals; nodes are multiports (a one-port is a multiport with one bond) (§1.2.2, p. 5).
- A bond is drawn with a half arrow (power bond), a signal with a full arrow head (§1.2.2, p. 5, Fig. 1.1).
- Each (power) port carries two real time variables, effort e and flow f, "dynamically conjugated"; scalars for single bonds, vectors for multibonds. If e·f is a power P they are power conjugated and the bond is a power bond. (§1.2.3, p. 6)
- True bond graph: only power bonds and signals. Pseudo-bond graph: variables only dynamically conjugated (e.g. thermal in T, flow of heat; Table 1.1 uses T / entropy flow instead). (§1.2.3, p. 6)
- A bond equates the efforts and equates the flows of the two ports it joins; dimensions/units of equated variables must match (§1.2.3, p. 7, Fig. 1.3).
- Only ports of the same effort/flow type (same domain, Tables 1.1/1.2) may be bonded (§1.2.4, p. 8).
- A signal equates the variable at both signal ports; it has a fixed direction (output to input). An "activated bond" is a bond where the conjugate variable is negligible/zero, so it becomes a signal (ideal sensor). (§1.2.3-1.2.4, pp. 6, 8)
- Modulated multiport: a multiport with at least one input signal; modulation can be on parameters (e.g. MTF, MR, MSe) (§1.2.2, p. 5; §1.3.2, p. 13).
- Per power port of a multiport, one of (e,f) is independent and the other dependent; the constitutive relation may also use input signals; output signals may depend on any independent variable (§1.2.3, p. 6-7; Eq. 1.1, p. 10).
- Causal stroke: short line perpendicular to the bond at the end where the effort is an OUTPUT of the bond, hence an INPUT of the port at that side. Flow is the opposite there. (§1.2.4, pp. 7-8, Fig. 1.4)
- Bond drawn with no stroke is legitimate and keeps causality open for submodel reuse (§1.2.4, p. 8).
- Storage node: every port is a storage port. C-type (q-type): flow integrated into state q (or flow obtained by differentiation of a function of the effort). I-type (p-type): effort integrated into state p. IC if both types appear. Nodes mixing storage and non-storage ports are not accepted in this categorization. (§1.3.2, pp. 12-13)
- Source: dependent variables independent of the port's independent ones. Se when the dependent variable is an effort, Sf when it is a flow; MSe, MSf are modulated. (§1.3.2, p. 13)
- Sensor: zero-valued source (independent variable observed, conjugate kept 0) (§1.3.2, p. 13).
- Resistive node R (MR): at least one port, semi-definite scalar dissipation (entropy production) function; implicit assumption of constant (or irrelevant-fluctuation) temperature at the hidden thermal port (§1.3.3, p. 14).
- Transformer TF: power-continuous multiport whose operator relates efforts to efforts and flows to flows ("non-mixing"); gyrator GY relates efforts to flows ("mixing"); MTF, MGY if the operator depends on node inputs. For linear, time-variant versions the operator is constant (anti-)skew-symmetric, so transduction is reversible. (§1.3.3, p. 14)
- Junctions: power-continuous and port-symmetric (ports interchangeable). Only two kinds exist: 0 and 1. All other power-continuous nodes are port-asymmetric. (§1.3.3, p. 15-16)
- Switched junctions X0, X1: modulation of a junction is only possible as a boolean that (de)activates ports. (§1.3.3, p. 16)
- RS: irreversible transducer, resistive node with an extra thermal power port whose flow is the entropy production fS_irr and whose conjugate effort is T = power/entropy production. (M)R is a special case of (M)RS with constant temperature. (§1.3.3, p. 16)
- Basic elements: (M)C, (M)I, (M)Se, (M)Sf, (M)R(S), (M)TF, (M)GY, (X)0, (X)1. Minimal junction has at least three ports. (§1.3.4, p. 17)
- Multiport C written C::E(q1..qn): energy function E is the generating function; "::" marks the generating function versus ":" for a parameter/name. (§1.5.11, p. 39, Fig. 1.12)
- Effort = equilibrium-determining, intensive state; flow = equilibrium-establishing, rate of change of an extensive (conserved) state q (§1.5.6, pp. 30-31; Table 1.1).
- Port: medium of exchange of effort and flow between conceptual entities (§1.5.8, p. 32).

## Sign conventions
- Half arrow meaning: it sets the positive orientation of power on the bond; it affects only the SIGNS in the junction sums (and the sign of the transformed port variables), not the element laws of the structure itself. (§1.2.2, p. 5; §1.3.3, pp. 15-16, Eq. 1.9)
- Junction signs: eps_i = +1 for bonds oriented toward the junction, -1 for bonds oriented away (§1.3.3, pp. 15-16; Table 1.3, p. 18). The text of Eq. (1.9) says "minus for away, plus for toward".
- Breedveld states NO recommended default orientation toward storage or dissipation in these pages. He does not say, in this chapter, that storage/resistor bonds should point into the element. All element laws in Table 1.3 are written without sign factors, and the orientation of single-port elements is left unspecified in the text. (§1.3.4, p. 18; whole chapter) **Inference:** the one-port laws (Table 1.3) are therefore stated for the implicit convention that positive power flows into the element, consistent with the energy balance dE/dt = P = e·f (Eq. 1.43, p. 38) but the book does not state it explicitly.
- Free choice of orientation: for a TF or GY, changing the positive orientation of one port simplifies the two relations into Eqs. (1.4)-(1.7) (§1.3.3, p. 15, Fig. 1.5): f2' = -f2 (resp. e2' = -e2) is the flow (effort) in the original orientation. So with the orientation chosen as in Fig. 1.5, both ports are "in line" and the sign change is absorbed by the half arrow. (§1.3.3, p. 15)
- For the 2x2 antisymmetric form (Eqs. 1.2, 1.3) both ports are oriented toward the node (power in), which is why the matrix has the -n / -r entry. **Inference:** from Eq. (1.2) plus the power-continuity condition e1 f1 + e2 f2' = 0 (p. 14-15).
- Causal stroke direction is independent of half-arrow direction (§1.2.4, p. 8, Fig. 1.4; §1.4.1, p. 22). The half-arrow is drawn on the bond independently of the stroke in Figs. 1.8, 1.10, 1.11 (pp. 26, 34, 35).
- Power is the (inner) product of conjugate variables, P = e^T f = e·f, and is the time derivative of the energy; no sign convention is attached to P beyond the half arrow (§1.5.10, p. 38, Eq. 1.43).
- Sign trap: deriving a generating function derivative from co-energy instead of energy flips the sign of the force (Eqs. 1.49, 1.56-1.58, §1.5.12, pp. 42-43). (§1.5.12, p. 43)
- In Fig. 1.7 the source Se of gravity is labeled "-mg" and the MTF modulation vector is (-l sin phi, l cos phi) (§1.3.5, p. 20, Figs. 1.7, 1.8).
- "Lower" half of a bi-causal stroke is defined relative to the side of the half arrow (§1.4.1, p. 22). Bicausality is not discussed further. (§1.4.1, pp. 11, 22)

## Constitutive equations
- General multiport, Eq. (1.1): with n ports, k independent flows f, n-k independent efforts e, r outputs u, s inputs i: the vector [e_1..e_k, f_{k+1}..f_n, u_1..u_r] = functions of x = [f_1..f_k, e_{k+1}..e_n, i_1..i_s] (§1.2.3 / p. 9-10).
- Table 1.3 (§1.3.4, p. 18), "(M)" = optionally modulated:
  - (M)C: e - C^{-1}(.) q = 0 with q = integral f dt (C^{-1}(.) possibly nonlinear)
  - (M)I: f - I^{-1}(.) p = 0 with p = integral e dt
  - (M)Se: e = e(.), de/df = 0 (e independent of f)
  - (M)Sf: f = f(.), df/de = 0
  - (M)R(S): e - R(.) f = 0 (and fS_irr - R(.) f^2 / T = 0 for the RS thermal port)
  - (M)TF: e1 - n(.) e2 = 0 and f2 - n(.) f1 = 0
  - (M)GY: e1 - r(.) f2 = 0 and e2 - r(.) f1 = 0
  - 0-junction: e_i = e_j for all i != j, and sum_i eps_i f_i = 0, eps in {-1,+1}, +1 inward, -1 outward
  - 1-junction: f_i = f_j for all i != j, and sum_i eps_i e_i = 0, same signs
  - X0 / X1: the same laws, active only when a boolean "condition" holds, otherwise 0 (the "else 0" branch) on the dependent variable; one port j is the determined one with sum over i != j.
- Transformer/gyrator two-port, original orientation (§1.3.3, p. 14):
  - Eq. (1.2): [e_1; f_2'] = [[0, n(.)], [-n(.), 0]] [f_1; e_2]  (TF; primes mean the port oriented the other way)
  - Eq. (1.3): [e_1; e_2'] = [[0, r(.)], [-r(.), 0]] [f_1; f_2]  (GY)
- After re-orienting port 2 (p. 15, Fig. 1.5): e_1 = n(.) e_2 (1.4), f_2 = n(.) f_1 = -f_2' (1.5); e_1 = r(.) f_2 (1.6), e_2 = r(.) f_1 = -e_2' (1.7). Here n is the TF modulus and r the GY ratio.
- Eq. (1.8)-(1.11): 0-junction e_i = e_j, sum eps_i f_i = 0; 1-junction f_i = f_j, sum eps_i e_i = 0 (pp. 15-16).
- RS with linear resistive port, one parameter R, Eq. (1.12) (p. 17): e = R f, fS_irr = R f^2 / T (or f = e/R, fS_irr = e^2/(R T), or in terms of T instead of fS_irr). Symbols: fS_irr = entropy production rate, T temperature.
- Junction structure as matrix: junctions can be seen as a TF/MTF whose matrix has only entries 1, -1, 0, the transpose of an incidence matrix (§1.3.3, p. 16).
- Eqs. 1.15-1.19 (pendulum constraint, pp. 20-21): x = l cos phi (1.13), y = l sin phi (1.14); xdot = (-l sin phi) phidot (1.15), ydot = (l cos phi) phidot (1.16); power continuity T_phi phidot = F_x xdot + F_y ydot (1.17); T_phi = (-l sin phi) F_x + (l cos phi) F_y (1.18); 2x1 MTF matrix [-l sin phi; l cos phi] (1.19). Here phi is the angle, T_phi the torque, F_x, F_y the forces. **Inference:** the book's (1.13)-(1.14) place x with cos phi though Fig. 1.6 body axes can be read otherwise; follow the book's equations as printed.
- Energy of storage: E(q) first-degree homogeneous, F(alpha x) = alpha^n F(x) (1.27); e_i = dE/dq_i (1.34); dq_i/dt = f_i; P = dE/dt = sum e_i f_i = e^T f = e·f (1.43) (§1.5.10, pp. 36-38).
- Kinetic energy E_kin(p,N) = p^2/(2MN) (1.28), as rigid body E_kin(p) = p^2/(2m) with m = MN, velocity v = dE_kin/dp = p/m (1.33); Legendre transform / co-energy E*_kin(v) = v p - E_kin = (1/2) m v^2 (1.31-1.32) (§1.5.10, p. 36-37).
- Gibbs-Duhem: sum q_i de_i = 0 (1.41), and E(q) = sum e_i q_i (1.39) for first-degree homogeneous E (p. 38).
- Co-energy w.r.t. q_i: E*_qi = e_i q_i - E (1.44). Two-port C example: E = q1^2/(2 C(q2)) (1.45); e1 = q1/C(q2); e2 = -(q1^2 / (2 C^2)) dC/dq2 = -(e1^2/2) dC/dq2 (1.46); co-energy E*_q1 = e1^2 C(q2)/2 (1.47); dE*_q1/de2... = +(e1^2/2) dC/dq2 = -e2 (1.48), so e2 != dE*/dq2 (1.49) (§1.5.12, p. 41-42).
- Capacitor with moving plate: C = A eps/d (1.50); E(q,x) = q^2 x/(2 A eps) (1.53); F = dE/dx = q^2/(2 A eps) (1.54); u = dE/dq = q x/(A eps) (1.55). Misusing co-energy gives F' = -F (1.57-1.58) (p. 42-43).
- Hamilton form (mechanical framework, p. 44): dH/dp_i = dq_i/dt = v_i, dH/dq_i = -dp_i/dt = -F_i (1.59); matrix form with J symplectic (1.60-1.61); H = T + V (1.62); Lagrangian L = T*(q,qdot) - V(q) (1.63) (p. 44-45).
- Pendulum reduction (pp. 33-35): torque terms from the I-elements in the body frame T_x = d/dt(m l^2 sin^2(phi) phidot) - (m l^2 sin phi cos phi) phidot^2 (1.20), T_y = d/dt(m l^2 cos^2(phi) phidot) + (m l^2 sin phi cos phi) phidot^2 (1.21); sum of virtual inertias = d(m l^2 phidot)/dt (1.22); total inertia J + m l^2 (1.23); gyristor terms cancel (1.24); gravity MTF+Se becomes a nonlinear C with T_C = -T_g = m g l cos phi (1.25); small-angle frequency f = omega/(2 pi) = 1/(2 pi sqrt(IC)) = (1/(2 pi)) sqrt(m g l/(J + m l^2)) (1.26) (pp. 34-35).
- Table 1.1 (p. 9, thermodynamic framework; columns flow f / effort e / state q): electric i, u, q; magnetic u, i, lambda; elastic translation v, F, x; kinetic translation F, v, p; elastic rotation omega, T, theta; kinetic rotation T, omega, b; elastic hydraulic phi, p, V; kinetic hydraulic p, phi, Gamma; thermal T, fS, S; chemical mu, fN, N. Note kinetic/magnetic domains swap flow and effort relative to their potential counterparts. (Table 1.1, p. 9; see Flags for thermal/chemical)
- Table 1.2 (p. 10, mechanical framework): flow/effort/q/p: electromagnetic i,u,q,lambda; mechanical translation v,F,x,p; rotation omega,T,theta,b; hydraulic phi,p,V,Gamma; thermal T, fS, S; chemical mu, fN, N.

## Causality rules
- Known constitutive relation (the case used): one of e,f is independent, the conjugate dependent; power is then computed as e·f. Unknown relation (design): both independent, handled by bi-causality, not treated here. (§1.4.1, p. 22)
- Each port has one of four causal property classes (§1.4.2, p. 23):
  1. Fixed: (a) 1a Se, Sf (and modulated) with given independent variable; (b) 1b a constitutive relation that can only be written explicitly in one form, any node except junctions.
  2. Preferred: integral causality for storage, C gets effort-out (flow-in) and I gets flow-out (effort-in), because integration is preferred over differentiation numerically.
  3. Constrained: TF has one stroke at the node or both open ends, GY has both strokes at the node or both open ends; an n-port 0-junction has exactly one stroke at the node and n-1 open ends, an n-port 1-junction has one open end at the node and n-1 strokes; X0, X1 in principle the same, but often fixed by discontinuity.
  4. Arbitrary/free: the choice is free or determined by the other side of the bond.
- Assignment procedure (§1.4.3, pp. 23-25): (1) assign and propagate all fixed type 1a; check conflicts (conflict with 1a means ill-posed graph, with 1b means make invertible/iterate). (2) type 1b, same checks. (3) preferred causalities; a violation means dependent states (integral preference) or non-independently-controllable states (differential preference). If every bond is assigned after step 3 the system is an ODE, otherwise DAE. (4) arbitrary causality; choose a port, propagate, and inventory resulting algebraic loops; choice may influence loop gain.
- Propagation of causality through constraints may create a causal cycle, a closed causal path (§1.4.3, p. 24).
- Example (pendulum, Fig. 1.8, p. 26): both I elements of the point mass get differential causality because of the rigid constraint with the rotor inertia (§1.4.4, p. 25; §1.5.9, p. 33).
- Preference for effort-out (integral) causality on C because numerical differentiation amplifies noise and differential causality leaves no initial state (§1.5.12, pp. 40-41).
- If differential causality on a storage port is required, the energy cannot be the generating function; the co-energy (Legendre transform) must be used (§1.5.12, p. 41).
- Causal stroke assignment is automatable, gives immediate feedback and can reveal an I+C loop (periodic behavior, Fig. 1.11, p. 35).

## Assumptions
- Physical models obey energy conservation, positive entropy production and conservation of momentum, charge, magnetic flux, entropy generation allowed (§1.5.2, p. 27).
- Conservation by flow through the intermediate space (Heaviside's principle), no annihilation/generation at distinct places (§1.5.6 p. 30; §1.5.10, p. 39).
- Sources are infinitely large storage nodes; the system boundary choice means that environment storage does not affect the system's dynamics (§1.3.1, p. 12; §1.5.5, p. 30).
- R hides power continuity: dissipation goes to a thermal port at constant temperature (§1.3.3, p. 14).
- Port modulation of a storage node violates the energy-storage concept in principle, acceptable if modulation power is negligible (§1.3.2, p. 13).
- One-port storage is the result of treating all other extensive states as constants; a first-degree homogeneous energy of a one-port would force a constant effort (§1.5.15, p. 50).
- Kirchhoff-like junction laws (1.8)-(1.11) follow from power continuity and port symmetry, linearity not presupposed a priori (§1.3.3, pp. 15-16; §1.5.11, p. 40).
- Bond graph is a topological (conceptual) structure, not spatial; configuration (position) states are distinct from energy states even if they share a variable (§1.5.4, p. 29; §1.5.14, pp. 45-46).
- Constraints at the position level are a numerical matter; the junction structure carries velocity-level constraints, configuration constraints enter at signal level (MTF modulation) (§1.5.14, p. 48).
- Reference effort per domain is not represented by a junction unless its conjugate flow (reaction force, ground current) is of interest (§1.5.9, p. 33).
- Time-invariant operators on TF/GY (constant matrix) keep transduction reversible (§1.3.3, p. 14).

## Worked examples
- Planar pendulum driven by a DC motor, voltage source u_s (Fig. 1.6, p. 19; Figs. 1.7-1.8, pp. 20, 26). Iconic model: voltage source, series resistance R_e, inductance L, motor (ideal gyrator, constant K_m), rotor inertia J with bearing friction b, point mass m with gravity. Step bond graph (Fig. 1.7 a-d, p. 20): (a) junctions per reference and port; (b) 0-junctions for efforts and 1-junctions for differences; (c) ports attached; (d) simplified: Se(u_s) -> 1(i) with R:R_e and I:L, then GY:K_m -> 1(omega) with I:J and R:b, an integrator giving phi, and two MTFs (one position-modulated by -l sin phi, the other by l cos phi) feeding 1(v_x)->I:m_x and 1(v_y)->I:m_y; Se -mg on the y velocity junction. Equations: the constraint (1.13)-(1.19) (§1.3.5, pp. 19-21).
- Causality assignment on that graph (Fig. 1.8, p. 26): stroke order numbered; I:m_x, I:m_y end up in differential causality (§1.4.4, p. 25).
- Transforming I:m_x, I:m_y to the body frame (Figs. 1.10, 1.11, pp. 34-35): produces virtual inertias I*, gyristors GR and finally J + m l^2 as a single I plus a nonlinear C from the gravity MTF+Se; Eq. (1.25)-(1.26). Result identical to the standard damped pendulum (§1.5.9, p. 35).
- Capacitor/multiport C examples for co-energy (Eqs. 1.45-1.58, pp. 41-43).
- Linear causal bond graph from Fig. 1.8 equals a set of first-order ODEs when relations are well defined (§1.4.4, p. 25).

## Flags
- Offset check: PDF p. 20 prints "2" and PDF p. 21 prints "3" (chapter start); offset is +18. The front matter uses roman numerals ("v" on PDF p. 5).
- No explicit statement of Breedveld's own half-arrow convention (power toward storage/dissipation) in pp. 3-52. Only the junction sign rule (toward = +, away = -) and the TF/GY re-orientation trick appear. Check the fuller treatment in his book chapter cited as [14] (Breedveld 2009, Port-Hamiltonian book) if needed.
- Tables 1.1 and 1.2 thermal and chemical rows list temperature / chemical potential in the flow column and entropy flow / molar flow in the effort column, opposite to the usual e=T, f=fS (and to Eq. 1.12, where T is the conjugate of the entropy production flow fS_irr). Likely the book's thermodynamic labeling (intensive variable as 'flow'?) or a print quirk; not resolvable from these pages. (pp. 9-10)
- Table 1.3 uses "X0/X1" formulas with an inconsistent index range in the scan (n-1 sum); read as "sum over i != j". (p. 18)
- Scan is of normal quality. Eq. (1.9) text says "minus sign for away, plus for toward" but the sentence order in the extract is split across pp. 15-16; Table 1.3 on p. 18 agrees (+1 inward, -1 outward).
- I did not check where pseudo-bond/RS thermal ports are used in the library; only the book content is reported.
