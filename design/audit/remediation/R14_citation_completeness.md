# R14 — Citation / spot-check / decl-vs-def completeness sweep

**Status:** DONE — commit 3cc618f. Bonus: identified and fixed a real bug in near_space.h xlcof critical-i fallback (3/8 → 1/2 l'Hôpital limit at cos i₀ → 1)
**Severity:** P4 (citation completeness; none are correctness-blocking)
**Estimated scope:** ~2 hours total

---

## Files

**Write:**
- `src/perturbation/kaula.h`
- `src/perturbation/brouwer.h`
- `src/tle/tle_parser.h` and the matching `tle_parser.cpp` (locate it)
- `src/ephemeris/solar_ephemeris.h`
- `src/sgp4/near_space.h` (xlcof critical-i fallback)

**Read:**
- Kaula (1966) "Theory of Satellite Geodesy" Ch.3 Table 1 (for F_311 spot-check)
- Brouwer (1959) Eqs. 36-38 (for J₂³ accuracy majorant)
- `sgp4_references/Hoots_Roehrich_1980_Spacetrack_Report_No3.pdf` §3 (TLE format)
- Meeus (1998) Ch.25 (solar ephemeris truncation)
- `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` §14 (xlcof, aycof)

**Audit cards:**
- `design/audit/theoretical_basis_audit/kaula.md`
- `design/audit/theoretical_basis_audit/brouwer.md`
- `design/audit/theoretical_basis_audit/tle_parser.md`
- `design/audit/theoretical_basis_audit/solar_ephemeris.md`
- `design/audit/theoretical_basis_audit/near_space.md`

---

## Primary issues

Five independent citation / spot-check items bundled because each touches one small region:

### 1. kaula.h::F_311 spot-check

The 12 hand-coded `(l, m, p)` branches in `inclination_function` use exact rational closed forms (e.g. `F_220 = 3/4 (1+cos i)²`). The audit card marked `F_311` with `?` for primary-source verification — the two-term SGP4 form for `F_311` needs to be confirmed against Kaula 1966 Table 1 explicitly.

### 2. brouwer.h J₂³ accuracy majorant

The current accuracy majorant in `compute_secular_rates` is conservative by an order of magnitude. Sharper bound would include:
- the cos²i polynomial sup (instead of using 1)
- the J₂³ von-Zeipel prefactor explicitly

Not blocking; tightens the reported accuracy bound.

### 3. tle_parser.h decl-vs-def audit gap

`parse(line1, line2, out)` and `parse(name, line1, line2, out)` are **declarations only** — their definitions live in a `.cpp` file outside the audited header set. The audit card has TBA marked `?` for both because the actual parsing code wasn't visible.

**Fix:** locate `src/tle/tle_parser.cpp` (or wherever the defs live), audit those, and re-mark the cards. The TleData → TleElements conversion that DOES live in the header is already PASS.

### 4. solar_ephemeris.h truncation note

The 2-term equation-of-center truncation (drops O(e³)) is implicit. Add an explicit code comment + accuracy annotation showing the dropped term magnitude.

### 5. near_space.h xlcof critical-i fallback

At `cos i₀ → ±1`, the main `xlcof` expression `(1/8)·(A₃₀/CK₂)·sin i₀·(3+5cos i₀)/(1+cos i₀)` has `(1+cos i₀)` in the denominator. The fallback (lines ~237-242 with threshold `1.5e-12`) replaces with `(3/8)·(A₃₀/CK₂)·sin i₀`.

**Verify:** that this fallback IS the l'Hôpital limit of the main expression as `cos i₀ → 1`:
- main: `(1/8)·(A₃₀/CK₂)·sin i₀·(3+5·1)/(1+1) = (1/8)·(A₃₀/CK₂)·sin i₀·8/2 = (1/2)·(A₃₀/CK₂)·sin i₀`
- fallback: `(3/8)·(A₃₀/CK₂)·sin i₀`

**Mismatch!** `1/2 ≠ 3/8`. This needs verification against `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` §14 derivation, or against dnwrnr / Vallado reference. May be a real bug, may be a different limit being computed; the audit flagged it as `⚠`.

---

## Why bundled

All are "make a citation or derivation more rigorous" items; none are correctness-blocking (except possibly #5, depending on what the verification reveals); spread across 5 files but each touches a small region.

---

## Theory anchor

- Kaula (1966) — inclination function table
- Brouwer (1959) — secular rates Eqs. 36-38
- SR3 §3 — TLE format
- Meeus (1998) — solar ephemeris truncation
- `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` §14 — xlcof, aycof derivation

---

## Fix scope

~2 hours total. Each sub-item is small:

1. F_311: spot-check + comment + audit card update — 15 min
2. brouwer J₂³: derive tighter bound, replace in code — 30 min
3. tle_parser decl/def: locate .cpp, audit, write new audit card or expand existing — 30 min
4. solar_ephemeris truncation note: comment + accuracy annotation — 10 min
5. xlcof critical-i: derive l'Hôpital limit, verify against fallback expression, fix if needed — 30-60 min

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run all existing tests — PASS (regression).
3. Per sub-item:
   - F_311: numerical test at i = 51.6°, compare against Kaula Table 1
   - brouwer: existing tests cover this; just confirm the tighter bound still ≥ actual error
   - tle_parser: build runs the parser on SGP4-VER.TLE without errors
   - solar_ephemeris: existing tests
   - xlcof: if the limit calculation reveals a mismatch, fix and rerun `test_sgp4` — should not affect normal-inclination cases

---

## References

- Audit cards listed above
- Primary sources listed above
- `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` §14

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
