# Louca et al., "Energy-Based Bond Graph Model Reduction", in Borutzky (ed.), Bond Graph Modelling of Engineering Systems: Theory, Applications and Software Support, Springer 2011: Ch. 2
PDF↔book page offset: book p. N = PDF p. N+18 (verified: book p. 53 = PDF p. 71, book p. 72 = PDF p. 90). Chapter occupies book pp. 53-103 (p. 104 blank).

Relevance: low. This chapter is about model reduction, which BondGraphToolkit.jl does not plan. Only items that may serve as diagnostics or reference examples are kept.

## Topics covered (one line each, with § and book page)
- Introduction and notion of a "proper" model for an application (§2.1, pp. 53-54).
- Metrics: activity, activity index, relative activity, Energetic Contribution Index (ECI) (§2.2, pp. 55-61).
- MORA, order reduction by dropping low-activity elements; quarter-car example (§2.3.1, pp. 62-72).
- Decoupling/partitioning by "conditioning" low-RA bonds into modulated sources; two-mass example (§2.3.2, pp. 72-82).
- ECI/SVD-based global reduction; slider-crank example (§2.3.3, pp. 82-88).
- Case studies: tractor semi-trailer, hybrid truck, pitch-plane truck, engine, HMMWV (§2.4, pp. 89-95).
- Discussion, limits, junction inactivity, software remarks (§2.5, pp. 96-100).
- Skipped: reduction flowcharts (Figs. 2.2, 2.4, 2.13, 2.14, 2.18, 2.19), tables of ranked results for the case studies, error plots, reference list, conclusion (§2.6).

## Definitions
- Bond power is the product of generalized effort and flow, P = e*f; the chapter's metrics are computed from it, not from state variables (§2.2, p. 55).
- Activity of an element over window [tau1, tau1+T]: A = ∫|P(t)| dt = ∫|e(t) f(t)| dt, the L1 norm of power; units of energy; unlike net stored energy it grows monotonically even over a full oscillation, where net energy change is zero (§2.2.1, p. 56, Eq. 2.2; Fig. 2.1, p. 57).
- Total activity A_total = Σ_i A_i over all k energy elements; activity index AI_i = 100 A_i / A_total (§2.2.1.1, p. 57, Eqs. 2.4-2.5).
- Relative activity of bond i at a junction: RA_i = A_i / A_max,junction, with A_max the largest bond activity at that junction. Below a user threshold epsilon the bond is "negligible": on a 0-junction its flow f_i drops out of the flow sum, on a 1-junction its effort e_i drops out of the effort sum (§2.2.2, p. 58, Eq. 2.6).
- Bond energy e_i(t) = ∫_{t0}^{t} effort*flow dτ; energy vector e(t) = [e_1 ... e_n]^T; energy trajectory matrix E stacks e^T(t_k) row-wise; energy Gramian W = ∫ e e^T dt ≥ 0 (§2.2.3, p. 59, Eqs. 2.7-2.9, 2.11).
- ECI (relative, per bond): c~ = Σ_i σ_i² |v_i| with (σ_i², v_i) eigenpairs of W; c = c~ / max(c~) (§2.2.3, p. 60, Eqs. 2.12-2.13).
- "Junction inactivity" (mentioned only): a 0- or 1-junction is inactive when all attached bonds have negligible activity (§2.5, p. 100).
- Subgraph, subgraph loop, driving/driven subgraph, driving/driven partition (definitions in §2.3.2, p. 74).

## Sign conventions
- Half-arrow = direction of algebraically positive power flow; a modulated source that replaces a conditioned bond keeps the same half-arrow direction as the original bond (§2.3.2, p. 73).
- Initial spring displacements in the slider-crank example: positive sign denotes extension (§2.3.3.1, p. 84).
- No junction sign equations are stated (they are used only implicitly: junction sums with positive/negative terms). Other convention facts: none in this chapter.

## Constitutive equations
- Element power from constitutive laws: I: P = e_I f_I = e_I Φ_I(p); C: P = Φ_C(q) f_C; R: P = e_R f_R = Φ_R(f_R) f_R, with Φ possibly nonlinear (§2.2.1, p. 56, Eq. 2.3).
- System form assumed: ẋ_i = F(x_i, x_d, u), x_i(0)=x_i0; 0 = G(x_i, x_d). The text says this covers derivative causality and Lagrange multipliers (§2.2, p. 55, Eq. 2.1).
- Conditioned (M)TF with a locally inactive bond is equivalent to an MSe with value n*e_1a (TF) or MSe n*f_1a (GY, shown as MSf then GY giving an MSe) (§2.3.2, p. 73-74, Table 2.6; table graphic is low resolution, see Flags).
- Quarter car (Table 2.1, p. 65): spring force a1 q + a3 q^3 for suspension and tire; tire force e(q) = a1 q + a3 q^3 when q > 0 (contact), e(q) = 0 otherwise (liftoff) (§2.3.1.1, p. 65).

## Causality rules
- Locally inactive bonds must be causally weak: they must not supply the flow input to a 1-junction nor the effort input to a 0-junction. Causality may need reassignment (§2.3.2, p. 73).
- Conditioning replaces a bidirectional power bond by a one-way signal, so the elements downstream are "driven" and do not back-excite upstream ones (§2.3.2, p. 73).
- Parasitic/constraint branches should not be conditioned; causality assignment is used to steer reduction so that the reduced model stays in explicit ODE form; reduction can introduce dependencies between storages (derivative causality) (§2.5, p. 98).
- Neglecting lever inertia in the two-mass example gives an algebraic loop; conditioning the c_2 output breaks it, cutting computation time by 35% (§2.3.2.1, pp. 78, 80-81).

## Assumptions
- Metrics are trajectory dependent: they hold for given inputs, initial conditions, parameters and time window; there is no proof that error is monotone in eliminated activity, only empirical evidence (§2.5, pp. 98-99).
- Reduced models are realization preserving (junction structure and physical meaning kept) (Abstract, p. 53; §2.3.1, p. 64).
- Elimination of elements with ideal-physical-element interpretation is done by a systematic procedure (ref. [24]) (§2.3.1.1, p. 68).

## Worked examples
- Quarter-car (sprung/unsprung mass, parallel spring-damper suspension and tire; 6 elements I,C,R x2; 4 states): bond graph in Fig. 2.5, p. 65; parameters Ms=267 kg, Bs=700 N s/m, Mu=36.6 kg, Bt=200 N s/m, suspension a1=18872 N/m, a3=2e5 N/m^3, tire a1=193915 N/m, a3=2e8 N/m^3; curb input Zr 0.2 m high, 0.5 m long, cycloidal; Se sources M g and a Sf for road velocity Vr(t) (§2.3.1.1, pp. 64-66, Table 2.1). Activity ranking at 5 m/s: suspension stiffness 45.44%, sprung mass 23.65%, susp. damping 11.70%, unsprung mass 10.61%, tire stiffness 8.457%, tire damping 0.1481%; beta=95% removes only tire damping (Table 2.2, p. 67). At 1 m/s the unsprung mass and tire damping are removed (Table 2.3, p. 69).
- Two-mass lever system (m1, m2, spring-dampers k1c1, k2c2, k3c3, lever a/b, step force F1): state-space with A, B matrices and one algebraic constraint (c3/c2) F_c2 - F_c3 + (b/a) c3 v1 - c3 v2 = 0 (Eq. 2.20); partitioned driving/driven equations in Eq. 2.21 (§2.3.2.1, pp. 78-81, Figs. 2.20-2.22).
- Slider-crank with rotational spring-damper, mass-spring-damper load and friction: kinematics v_G, v_C and constraint v_Cy = L1 α̇ cos α + L2 β̇ cos β = 0 (Eqs. 2.22-2.23); bond graph with 31 numbered bonds in Fig. 2.25; parameters in Table 2.7 (§2.3.3.1, pp. 83-84).
- DC motor: none in this chapter.

## Ideas usable as library diagnostics/tests
- **Inference:** per-bond power `P_i(t) = e_i f_i` (using the half-arrow orientation) is trivially available from simulation outputs; a helper exposing it is the base for everything below (§2.2, p. 55).
- **Inference:** power-balance test: for each junction, Σ_i ± P_i should be zero for all t (follows from equal-effort/equal-flow plus the sum constraint, not stated as a test in the chapter); and for a closed system, Σ_elements P = 0 with dE_storage/dt = Σ P_R-sources etc. (§2.2, p. 55; §2.2.1, p. 56).
- Element activity A_i = ∫|e f| dt and activity index AI_i as a per-element energy-throughput diagnostic; cumulative AI sorted descending ranks elements by importance (§2.2.1, pp. 56-57; §2.3.1, p. 63, Eq. 2.18).
- Relative activity per junction bond: cheap check of which bonds are weakly coupled (§2.2.2, p. 58).
- Cheap sanity test from an example: for the C/I elements, ∫P dt over a full oscillation cycle is about 0 while activity is strictly positive (§2.2.1, p. 56; Fig. 2.1, p. 57).
- ECI recipe using SVD of S with S_kj = sqrt(Δt_k) E_kj; eigenvalues of W are the squared singular values (§2.2.3.1, p. 61, Eqs. 2.14-2.16; §2.3.3, p. 82).
- Reference models with published numbers (quarter-car activity tables, two-mass system) could be test fixtures: e.g. quarter-car activities in J at 5 m/s: 4139, 2155, 1066, 966.2, 770.7, 13.49 (Table 2.2, p. 67). **Inference:** reproducing them would require the same hand-derived equations and curb input, so use only as loose cross-checks.
- Accuracy metric: average relative error ε = 100 ∫|w - w_r| dt / ∫|w| dt for comparing a model against a reference (§2.3.1.1, p. 68, Eq. 2.19).

## Software remarks
- 20-Sim's element library has 1- and 0-junctions that compute the relative activity of connected bonds; AMESim implements activity analysis and reports results after simulation, with reduction done manually (§2.5, p. 96).
- Quarter-car time responses were generated in a bond graph environment (20-Sim, CAMP-G) from hand-derived equations transferred to MATLAB (§2.3.1.1, p. 66).
- Metrics can also be computed in non-bond-graph tools (ADAMS, MATLAB) (§2.5, p. 96).

## Flags (anything unreadable, ambiguous, or not verifiable)
- Tables 2.5 and 2.6 (p. 73-74) are low-resolution diagrams; the TF/GY conditioning equivalences above are read from the figure and the prose, treat with caution.
- Fig. 2.25 (p. 84) and Fig. 2.34 (p. 92, rotated) are too small to transcribe junction/MTF modulations reliably; not recorded.
- The chapter never states junction signs or a formal causality-assignment algorithm; it assumes standard SCAP/derivative causality knowledge.
- Eq. 2.20 matrix entries were read from a scan; verify against the page before using as a test fixture (p. 79).
