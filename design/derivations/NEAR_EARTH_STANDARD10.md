# Near-earth assembly — per-line Standard-10 coverage (AUDIT_BACKLOG A2)

Every boxed formula in the near-earth propagation *assembly* (not just the drag coefficients, which are
already code-matched) gets a verifier asserting `simplify(code − reference) == 0` against the canonical
SGP4 reference (Spacetrack Report #3 / Vallado SGP4 / dnwrnr / python-sgp4). This is roadmap items W1–W5;
the **W5 Definition of Done is that the "unverified" count below is 0.**

| Stage | Boxed formulas | Verifier | Status |
|---|---|---|---|
| Short-period J₂ periodics (`src/perturbation/short_period.h`) | `rk`, `uk`, `xinck`, `xnodek`, `rdotk`, `rfdotk` (6) | `verify_near_short_period.m` | ✅ **W1** (6/6 PASS) |
| Long-period / Lyddane recovery (`src/sgp4/near_space.h`) | `xlcof`/`aycof` coefficients; `axN`, `ayN`, `xll`, `xlt` | `verify_near_longperiod.m` | ✅ **W2** (5/5 PASS) |
| Kepler solve (`src/orbit/modified_kepler.h`) | Newton + Halley residual `f`, denominator `fp=−f'`, increments | `verify_near_kepler.m` | ✅ **W3** (6/6 PASS) |
| Secular advance + drag assembly (`src/orbit/secular_update.h`) | `tempa`/`tempe`/`templ`; `delm` η-cubic; `am`/`em`/`nm`; `nodecf·t²` | `verify_near_secular.m` | ✅ **W4** (11/11 PASS) |
| Osculating → TEME position (`src/orbit/state_from_elements.h`) | `ux/uy/uz`,`vx/vy/vz` orientation (+ orthonormality proof); `r·U`; velocity | `verify_near_osculating.m` | ✅ **W5** (9/9 PASS) |

**Unverified stages: 0 — all near-earth assembly stages are Standard-10 verified (P1 COMPLETE).**

Together with the already-code-matched drag coefficients and the deep-space (SDP4) subsystem, the entire
SGP4 path is now provable formula-by-formula: authentic-SGP4 fidelity is no longer merely end-to-end
(33/33) but per-line `simplify(code − reference) = 0`.

> W3 also corrected a sign-flipped Kepler equation in `modified_kepler.h`'s file docstring
> (it read `x − ayn·cos x + axN·sin x = U`; the code and Vallado solve `x + ayn·cos x − axN·sin x = U`).

## Notes

- `half_J2` in the code is `CK2 = J2/2`; the reference's `temp1 = 0.5·J2/p` equals the code's
  `sp2 = half_J2/p` exactly. Verifiers encode both sides from their own intermediate definitions so a
  coefficient drift on either side breaks the identity.
- These verifiers are *transcription guards*: the near-earth assembly implements the published SR3 closed
  forms directly (unlike the drag coefficients, which required multi-step derivation). The check proves the
  code reproduces the canonical equation and stays a regression guard against future edits.
- Run a single verifier: `octave-cli --no-gui --quiet design/derivations/verify_near_short_period.m`
  (exit 0 = PASS). The full set runs under `tools/run_acceptance.ps1` (gates W1–W5).
