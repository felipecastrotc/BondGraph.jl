# Ticket 02 report

Implemented and exported `Mass`, `Spring`, and `Damper` as mechanical shorthands for the role-named one-ports, including the v0.2 defaults. Ported the 1-DoF and Duffing Literate examples and golden regression models to `Inertia`, `Compliance`, and `Resistance`, with the source bond into the 1-junction and element bonds out. Added a short one-port Literate example and behavioral tests for shorthand trajectories, inertia energy, and two inertias on a 1-junction.

Tests: offline Julia 1.10.12 `Pkg.test()` passed, 43/43 assertions. The four mechanical golden cases passed without changing their CSVs. All three changed or new Literate files parsed; the 1-DoF example's named model helper compiled. The docs build was outside this ticket's scope.

Surprise: an effort-driven two-inertia junction compiles to a mass-matrix DAE that `Tsit5` cannot solve. The two-inertia simulation check uses a time-varying flow source and verifies both momenta against the closed form, without adding a solver dependency.

Open questions: none for ticket 02. A future DAE solver would allow an effort-driven two-inertia trajectory check.

Commit blocker: `git add` could not create `/Users/fctc/me/doutorado/dev/bg/bg_lib/.git/worktrees/bg_lib-t02/index.lock` (`Operation not permitted`), even though the worktree files are writable. The changes and this report remain uncommitted.

STATUS: BLOCKED
