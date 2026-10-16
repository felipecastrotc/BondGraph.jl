# Ticket 03: multi-junction models

- Ported the 2-DoF and four-mass N-DoF Literate examples and golden regressions to role-named elements and named free ports joined with `connect`.
- Drew the N-DoF coupling bonds `0 → sd`, initialized the second coupling compliance at `q = -1`, and flipped the coupling `q` golden columns in the regression mapping.
- Added a public API check for the closed-form RC trajectory, junction power balance, and flow balance across a free-port bond. Removed both cases from `PENDING`.

Test result: offline Julia 1.10.12 `Pkg.test()` passed, 43/43 tests. The first red run found that MTK's exposed free-port flow sign is opposite the bond-flow sign used in my initial power assertion; the corrected assertion passes. No changes to the core or golden CSVs.

Open questions: none for the implementation. The docs build remains with ticket 10 and was not run here.

Commit blocked: `git add` could not create `.git/worktrees/bg_lib-t03/index.lock` (`Operation not permitted`). A direct write probe to that Git directory failed too, so the changes remain unstaged in this worktree.

STATUS: BLOCKED
