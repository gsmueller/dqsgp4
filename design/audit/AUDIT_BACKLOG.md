# Standards-Compliance Audit Backlog

> Created 2026-06-03 (maintenance/accounting pass). Purpose: keep current information & code
> "keeping up" with the standards the project already set for itself (Standard 10, the charter→
> spec→audit→test traceability in `design/index.md`, the zero-magic-numbers policy, the error
> framework, the sealed-room/cleanroom derivation discipline). Each item is a GAP where current
> work has not yet been held to an established standard — not new scope.
>
> Grounding facts (verified this pass, not assumed): all **75 `verify_*.m`** pass; the SGP4
> **near-earth drag** coefficients are Standard-10 code-matched (`simplify(derived−code)=0`, 112
> drag-family checks); the **deep-space/SDP4** subsystem is **not** derived or audited to Standard
> 10 (9 failing sats); `design/index.md` coverage = specs 2/6, audit layer 3/5, tests 0/5;
> `design/documentation_gaps.md` lists opaque constants; REQ-EF-6 (model-truncation accuracy) is
> marked a gap in `audit/test_coverage.md`. Confidence is flagged per item.

## P1 — blocking / highest value

### A1. Deep-space (SDP4) held to Standard 10 — and root-cause the 9 failures
- **Checks:** every deep-space coefficient/term (`dscom`/`dsinit`/`dspace`/`dpper` analogues in
  `src/sgp4/deep_space.h`) satisfies `simplify(derived−code)=0` against the born-digital
  `python_sgp4_rhodes` reference, exactly as the near-earth drag coeffs were; and the 9 failing
  deep-space sats reach the 0.1 km deep-space tolerance.
- **Enforces:** Standard 10 (the standing code-match rule); it is the one SGP4 subsystem not yet
  at the bar the drag work set. **Scope:** `src/sgp4/deep_space.h`, `tests/test_sgp4/`. **Effort:** L.
- **Confidence: HIGH.** This is the documented next thread; the asymmetry (drag audited, deep-space
  not) is real. Pursue per the handoff's 3 root-cause classes (A lunar-solar / B resonance / C error-code).
- **Progress (2026-06-03, class C RESOLVED → 24/33 → 25/33):** sat 33334 was the sole t=0 error case.
  The reference returns error 3 (perturbed ecc out of range) with `ep=-122.217193`; our `deep_space.h`
  *detected* it but silently returned `(0,0,0)` instead of an error, so the harness scored a 45005 km
  miss vs the reference's stale buffer. Fix: NaN error-state (matches the reference r=(NaN,NaN,NaN)),
  reference-exact guard `e<0||e>1` (propagation.py:1837), stale-buffer detection in `tcppver_parser.h`
  (mirrors python-sgp4 `"(Use previous data line)"`), and expected-error scoring in `main.cpp`.
  **Standard-10 win for class A:** our DPPER perturbed `e` matched the reference's `-122.217193` exactly,
  so the deep-space DPPER *eccentricity* path is correct. Verifier:
  `tests/test_sgp4/verify_deep_space_error_codes.py` (5/5, locks the canonical `[1,1,6,6,4,3,6]`).
- **Progress (2026-06-03, classes A+B RESOLVED → 25/33 → 33/33 — FULL SUITE PASSES):** the handoff's
  "resonance integrator" framing for class B was **refuted** — all four (28623/23599/11801/16925) are
  `irez=0` (non-resonant), high-e, high-**drag** deep-space orbits. The real bug: `propagate_deep_space`
  **omitted the near-earth secular drag** (`tempa`/`tempe`/`templ`, `nodecf·t²`, and non-simple
  `delomg`/`delm`) that the reference applies to ALL satellites (propagation.py:1713-1794) — our
  deep-space path was effectively drag-free, which is why the only passing high-e deep-space sat (20413)
  has `bstar=0`. Zeroing drag in the reference reproduced our error to every digit. Fix: apply the
  near-earth drag in the deep-space path in the reference order (drag→DSPACE→`am=(xke/nm)^⅔·tempa²`,
  `em−=tempe`, `M+=n0·templ`), threading the always-computed `NearSpaceInit` drag coeffs into
  `propagate_deep_space`. This **also fixed class A** (22674/23177/21897/23333 were the same missing
  drag — dominant over the lunar-solar drift the handoff hypothesized). Verifier:
  `tests/test_sgp4/verify_deep_space_drag.py` (13/13). **test_sgp4 = 33/33, 623/623.**
- **Residual ROOT-CAUSED + FIXED (2026-06-03c, commit pending):** the class-A residual was the
  **simple-vs-non-simple drag (`isimp`) flag**. `sgp4init` forces `isimp=1` (SIMPLE drag) for EVERY
  deep-space sat (period ≥ 225 min, propagation.py:1496-1499), independent of perigee — but
  `propagate_deep_space` gated the non-simple drag on `ns.use_simple_model`, a *perigee < 220 km* test.
  So deep-space sats with perigee > 220 km (23177 @ 349 km, 22674/21897/23333) **wrongly received the
  non-simple drag**, whose `bstar·C5·sin(M)` term is periodic → an oscillating eccentricity error.
  Localized by element-level instrumentation (post-DPPER mean elements vs `python_sgp4_rhodes`): the diff
  was entirely in `em` (≤1.27e-7, oscillating), inc/node/argp matched to ~1e-13. Fix: deep-space path
  uses the **simple** drag model only. Result: **23177 6.1e-3 → 7.1e-9 km; 22674 1.9e-2 → 3.9e-6 km;
  21897 4.1e-4 → 1.1e-6 km** (`em` now matches the reference to machine precision). Verifier
  `verify_deep_space_drag.py` extended (15/15) to lock the isimp invariant (18 deep-space sats have
  perigee>220 but isimp=1). Earlier ruled-out (still true): secular rates match to ~2e-16, DPPER
  coeffs/formula/constants correct, the pre-DPPER mod-2π reduction is inert.
- **Resonance residual ROOT-CAUSED + FIXED (commit `4496fda`) — the WGS72 earth-rotation rate.** The
  `irez=2` half-day-resonant sats (22674/21897) had a ~1e-6 km residual. Element-level instrumentation
  localized it: post-resonance `xn` matched to 1.5e-15, but the mean longitude `xll` was off by 6.0e-11 rad.
  Cause: the sgp4_standard WGS72 `omega` was truncated to `7.292115e-5` (WGS84's value); the deep-space
  resonance `theta_dot = omega·60` must equal the canonical SGP4 `rptim = 4.37526908801129966e-3` rad/min
  (propagation.py:672). The full-precision WGS72 `omega = 7.2921151467e-5` matches `rptim` to 8.7e-15
  (vs an 8.8e-11 gap). `omega` is used ONLY for the resonance `theta_dot`, so no other sat is affected.
  Result: **22674 3.9e-6 → 1.5e-7 km; 21897 1.1e-6 → 2.3e-8 km.** `verify_deep_space_drag.py` 17/17.
- **23333 is the double-precision FLOOR (not a bug; ~8.8e-6 km, far sub-tolerance).** Extreme orbit
  (e=0.973 mean / 0.99 perturbed, a=37.9 ER). Mean elements match the reference to ≤7e-12; the Kepler
  converges identically (both the reference Newton-10-cap+limiter and a full solve give the same E in
  8 iters — the limiter never fires — so our Halley solver also matches). The 8.8 mm is the ~7e-12
  rounding differences in the DPPER/assembly amplified by the extreme-e position sensitivity (~1.6e6 km/rad).
  Reducible only by bit-exact operation-order matching; not worth it. **All deep-space residuals now ≤1.5e-7 km
  except this one extreme-e floor.**

### A2. Near-earth SGP4 propagation assembly — Standard-10 per-line, not just coefficients
- **Checks:** the drag *coefficients* are code-matched, but the near-earth *propagation* path
  (short-period, long-period/Lyddane, Kepler solver, osculating→position assembly in
  `near_space.h`/`secular_update.h`/`sgp4_propagator.h`) has no `simplify(derived−code)=0` verifier
  per line — only end-to-end output validation (≈1e-9 km on the passing sats, which is strong).
- **Enforces:** Standard 10. **Scope:** the non-coefficient near-earth code. **Effort:** M–L.
- **Confidence: MEDIUM.** The output validation already gives high empirical assurance; this closes
  the *formal* per-step gap. Lower urgency than A1.

## P2 — completeness of the established framework

### A3. Finish the specification layer (REQ-SY / REQ-CP / REQ-PR / REQ-IN)
- **Checks:** `design/index.md` marks `specifications/{system,constants_provider,propagator,integrator}.md`
  as *forthcoming* (2/6 written). The DQ propagator + constants provider are live code with only
  implicit requirements. **Enforces:** the index.md traceability invariant (every REQ cites an OBJ/CON;
  every test cites an AUD which cites a REQ). **Scope:** `design/specifications/`. **Effort:** M.
- **Confidence: HIGH** (index.md states it outright).

### A4. Finish the audit layer + stand up the test layer (AUD-DOC, AUD-TC; tests 0/5)
- **Checks:** `audit/documentation.md` (AUD-DOC) and `audit/test_coverage.md` (AUD-TC) are forthcoming;
  the index.md test layer is 0/5 implemented. Without AUD-TC the AUD↔test mapping can't be enforced.
- **Enforces:** index.md layers 3–4. **Scope:** `design/audit/`, `tests/audit/`. **Effort:** S–M.
- **Confidence: HIGH.**

### A5. Force-model accuracy propagation (REQ-EF-6 / AUD-EF-6) — close the marked gap
- **Checks:** every force lambda (`gravity_central`, `gravity_zonal`/J2, `drag`, future third-body)
  adds its model-truncation residual into `state.errors.accuracy`. `test_coverage.md` marks this a gap;
  only `drag.h` partially does it. **Enforces:** REQ-EF-6, OBJ-3 (three-error tracking). **Scope:**
  `src/forces/`. **Effort:** S–M. **Confidence: MEDIUM-HIGH** (the gap is explicitly recorded).

## P3 — quality / drift / lower urgency

### A6. Zero-magic-numbers compliance sweep
- **Checks:** every numeric literal in `src/` traces to a `design/zero_magic_numbers_policy.md` category;
  `design/documentation_gaps.md` flags ~dozens of opaque constants (Hansen/Kaula/ephemeris/resonance).
- **Enforces:** the zero-magic policy, AUD-CC-15. **Scope:** `src/`. **Effort:** M. **Confidence: MEDIUM**
  (real, but much of it is the BH61-textbook constant-provenance effort, partly tracked elsewhere).

### A7. Verifier + dimensional-audit completeness for non-drag boxed formulas
- **Checks:** every boxed formula in the live derivations has both a Standard-10 check and a
  dimensional-audit check (the drag work has both; Brouwer/Kaula/element-recovery may not).
- **Enforces:** Standard 10 + the dimensional-audit rule (`feedback_dimensional_audit`). **Scope:**
  `design/derivations/` ↔ `verify_*.m`. **Effort:** M–L. **Confidence: MEDIUM.**

### A8. DQ-algebra symbolic verification (Standard 10 for REQ-DQ ops)
- **Checks:** the dual-quaternion algebra laws (REQ-DQ-1..18) have symbolic (not just unit-test)
  verifiers — pose composition, screw exp/log round-trip, duality. `test_dual_quaternion` is 72/72
  but unit-level. **Enforces:** Standard 10 applied to DQ. **Scope:** `src/math/dual_quaternion.h`.
  **Effort:** M. **Confidence: MEDIUM.**

### A9. Code↔documentation drift re-validation
- **Checks:** `design/code_to_documentation_xref.md` still covers every constant/formula in code after
  the recent additions (DQ integrator, new force lambdas). **Enforces:** OBJ-7, AUD-CC-3. **Scope:**
  the xref doc vs `src/`. **Effort:** M. **Confidence: LOW-MEDIUM.**

## Notes on provenance
This list was seeded by an Explore agent over the governance docs, then **vetted** against
verified facts (the agent had repeated the stale "89 verifiers" figure — corrected to 75 files /
112 drag checks — and named a few files I could not confirm exist, which were dropped). Treat
effort/confidence as estimates. A1 and A3/A4 are the highest-leverage: A1 is the live technical
thread; A3/A4 unblock the project's own traceability and test story.
