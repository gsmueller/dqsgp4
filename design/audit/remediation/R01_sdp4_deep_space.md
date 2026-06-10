# R01 — SDP4 deep-space implementation

**Status:** PARTIAL — 14/24 deep-space PASS (was 0/24); 10 remaining FAILs are high-B* drag-dominated or Lyddane edges
**Severity:** P0 (blocking C-FAIL)
**Estimated scope:** ~1 work-day (largest of 14)
**Origin:** AUD-TBA findings; coupled C-FAIL #1 + #2

---

## Files

**Write:**
- `src/sgp4/deep_space.h`
- `src/perturbation/resonance.h`

**Read:**
- `sgp4_references/dnwrnr_sgp4/libsgp4/SGP4.cc` (reference C++ implementation — has DPPER, DSPER, d-coefficients)
- `sgp4_references/dnwrnr_sgp4/libsgp4/SGP4.h`
- `sgp4_references/vallado_celestrak/software/cpp/SGP4/SGP4/SGP4.cpp` (Vallado reference)
- `sgp4_references/Hoots_Roehrich_1980_Spacetrack_Report_No3.pdf` (primary source — §6 SDP4)
- `tests/test_sgp4/main.cpp` (validation harness)
- `sgp4_references/aholinch_sgp4/data/tcppver.out` (24 deep-space reference points)

**Audit cards:**
- `design/audit/theoretical_basis_audit/deep_space.md`
- `design/audit/theoretical_basis_audit/resonance.md`

---

## Primary issue (C-FAIL #1 + #2 coupled)

### deep_space.h

`propagate_deep_space` cites SR3 §6 but omits three pieces:

1. **DPPER** — Sun + Moon long-period periodic corrections to e, i, M, ω, Ω
2. **DSPER** — 24h and 12h tesseral resonance integration (which consumes `ds.resonance`)
3. **Time-varying sun/moon angles** in the secular rates

The function currently reduces to "near-space with constant solar/lunar offsets" — not enough to match the SR3 reference output for any of the 24 deep-space test cases.

### resonance.h

`step_resonance` is a stub:
- `(void)earth_rate` discards the input
- Uses linear extrapolation at constant `n` instead of SR3 §6 tesseral acceleration
- At initialize, 10 d-coefficients (d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433) are silently zero-errored (TODO at `resonance.h` line ~134)

---

## Why coupled

`deep_space.h::propagate_deep_space` calls into `resonance.h::step_resonance`. They share:
- The SDP4 silo (SR3 §6)
- The same test vectors: 24/24 deep-space `tcppver.out` cases
- The same reference implementation pattern (DPPER + DSPER + secular update with time-varying angles)

Splitting would create a serialization point with rework risk.

---

## Theory anchor

Hoots-Roehrich 1980 Spacetrack Report #3 §6 (SDP4):
- DPPER: long-period Sun/Moon periodic corrections (closed-form)
- DSPER: tesseral resonance integration (leapfrog or Euler at 720 min)
- Time-varying angles: GMST + Sun longitude + Moon longitude rates

dnwrnr libsgp4 (`SGP4.cc`) implements all three; output matches `tcppver.out` for all 24 deep-space cases.

---

## Fix scope

Approach:

1. **Port DPPER** from `dnwrnr_sgp4/libsgp4/SGP4.cc` (Sun + Moon long-period periodic corrections). Apply to e, i, M, ω, Ω after secular advance, before short-period.

2. **Populate the 10 d-coefficients** in `initialize_resonance` per SR3 §6 formulas:
   - d2201, d2211 — 24-hour geosynchronous resonance
   - d3210, d3222, d4410, d4422 — same
   - d5220, d5232, d5421, d5433 — 12-hour Molniya resonance
   - dnwrnr provides exact formulas

3. **Implement integration** in `step_resonance`:
   - Choose leapfrog (SR3 specifies 720-min step) or Euler
   - Add REQ-EF-5 residual to precision per `kepler.h:79` pattern
   - Remove the `(void)earth_rate` discard; use the rate

4. **Wire `propagate_deep_space`** to call:
   - Secular advance (with time-varying angles)
   - Tesseral resonance integration (DSPER)
   - DPPER long-period periodic corrections
   - Short-period (already correct)

---

## Verification

Build + `test_sgp4.exe` against `sgp4_references/aholinch_sgp4/data/tcppver.out`:

- **Expected: 0/24 → 24/24 deep-space PASS**
- Regression: 8/8 near-earth cases must still PASS at ~7e-9 km max position error

Per-satellite expected position-error thresholds for deep-space (per dnwrnr's match rate):
- sat 04632 (period ≈ 1198 min): < 0.1 km
- sat 23333 (long-arc, 17600 km error currently): < 1 km
- All 24: < 0.1 km maximum position error

---

## References

- AUD-TBA audit cards: `design/audit/theoretical_basis_audit/deep_space.md`, `resonance.md`
- Consolidated summary: `design/audit/AUD_TBA_results.md` §"The three C-FAILs"
- SR3 derivation transcript: `sgp4_references/hoots_roehrich_1980/hoots_roehrich_1980_math_derivation.md`

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
