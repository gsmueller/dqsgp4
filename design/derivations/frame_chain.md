# The Ecliptic-of-Date → GCRS Frame Chain

*Theory note preceding the precession/obliquity rotations in `astronomy/frames.h` and the
`body_position_gcrs` consumer (L3 ephemeris frame chain). Per the theory-first mandate the
chain's rotations are derived here so the functions fall out; per the no-perceived-fidelity
mandate the independent oracle is fixed before code. Sources, all born-digital / in-repo:
Capitaine, Wallace & Chapront (2003) A&A 412, 567 (IAU 2006 precession); IERS Conventions
(2010) §5.4–5.6; ERFA `pfw06.c`/`pmat06.c`/`obl06.c`/`ecm06.c` (the SOFA reference port);
coefficients pinned in `born_digital_sources.md` §A (VERIFIED bit-exact vs erfa.pfw06/pmat06).*

## 1. The problem

`sun_meeus` / `moon_meeus` produce a body's geometric position as ecliptic longitude/latitude
λ, β and distance Δ, in the **mean ecliptic and equinox OF DATE** (Meeus §25/§47). Any
consumer that works in the inertial celestial frame — the SR3 lunisolar perturbations, an
observer reduction, a comparison against a JPL ephemeris — needs that position in **GCRS**,
the Geocentric Celestial Reference System (the geocentric realization of the ICRS, ≈ J2000
mean equatorial). The transform is a fixed chain of rotations parameterized by time t (TT
Julian centuries from J2000.0):

    ecliptic-of-date  --Rx(−ε_A)-->  equatorial-mean-of-date  --P(t)ᵀ-->  GCRS.

Two rotations: the obliquity rotation (ecliptic→equatorial) and the inverse bias-precession.

## 2. Ecliptic spherical → Cartesian

    r_ecl = Δ · ( cos β cos λ,  cos β sin λ,  sin β ).                          (1)

A unit-aware Cartesian triple in the ecliptic-of-date frame; the rotations below are pure
orthogonal Matrix3 (L2 `rot_x`/`rot_z`), so Δ's units pass through unchanged.

## 3. Obliquity rotation — ecliptic → equatorial of date

The ecliptic and the equator of date share the x-axis (the equinox γ); the equator is tilted
from the ecliptic by the **mean obliquity** ε_A(t). Rotating a vector from ecliptic to
equatorial coordinates is a rotation about x by −ε_A:

    obliquity_rotation(ε_A) ≡ Rx(−ε_A),     r_equ = Rx(−ε_A) · r_ecl.          (2)

ε_A is the IAU 2006 mean obliquity (IERS Conv. 2010 Eq. 5.40) — exactly what
`astronomy/obliquity.h::obliquity_iau2006` already generates and which is currently an
ORPHAN (no consumer). This chain wires it. (Sign check: ERFA's `ecm06` builds the GCRS→
ecliptic map as Rx(ε_A)·P, so equatorial→ecliptic is Rx(+ε_A) and its inverse ecliptic→
equatorial is Rx(−ε_A) — born_digital §A records ecliptic→equatorial = Rx(−ε) verified exact
at t=0 vs astropy.)

## 4. Bias-precession — equatorial-of-date → GCRS

The mean equator/equinox of date are carried back to GCRS by the inverse IAU 2006
**bias-precession** matrix P(t) (GCRS → mean equatorial of date; the "06" matrix folds in the
frame bias ICRS↔J2000, IERS Conv. 2010 §5.4). P is parameterized by the Fukushima-Williams
four-angle set (γ̄, φ̄, ψ̄, ε_A) — Capitaine+ 2003 / `pfw06` — via the composition (§A)

    P(t) = fw2m(γ̄, φ̄, ψ̄, ε_A) = R1(−ε_A) · R3(−ψ̄) · R1(φ̄) · R3(γ̄),          (3)

with R1 = `rot_x`, R3 = `rot_z`. The four FW angles are secular polynomials in t (arcsec,
×π/648000 → rad), pinned born-digital in §A:

    γ̄(t) = −0.052928   + 10.556378·t   + 0.4932044·t²  − …
    φ̄(t) = 84381.412819 − 46.811016·t   + 0.0511268·t²  + …
    ψ̄(t) = −0.041775    + 5038.481484·t + 1.5584175·t²  − …
    ε_A(t) = obl06(t)   (the same mean obliquity as §3 — obliquity.h)

Each coefficient is a finite-digit IAU fit → `model_coefficient` (digits → accuracy, binary
storage → T-scaling precision), evaluated through a `TrackedPolynomial` exactly as
`obliquity.h` treats Eq. 5.40. So the precession matrix carries the three-error budget and
its precision tightens with a wider T. Equatorial-of-date → GCRS is Pᵀ.

## 5. The composite

    r_gcrs = Pᵀ · Rx(−ε_A) · r_ecl  ≡  M(t) · r_ecl,    M(t) = P(t)ᵀ · Rx(−ε_A).   (4)

M is precisely the transpose of ERFA's ecliptic-to-ICRS matrix `ecm06` (since ecm06 =
Rx(ε_A)·P ⇒ ecm06ᵀ = Pᵀ·Rx(−ε_A) = M), which gives a single element-wise oracle for the whole
chain.

**Fidelity — mean-consistent, no nutation term.** Meeus §47/§25 give the body referred to the
**mean** equinox/ecliptic of date, and ecm06 = Rx(ε_A)·P is exactly the **mean**-ecliptic↔GCRS
rotation (mean obliquity + bias-precession, no nutation). So the chain is mean-CONSISTENT:
rotating a mean-ecliptic position by the mean rotation Pᵀ·Rx(−ε_A) is the EXACT coordinate
transform of one physical point, recovering the body's GCRS position to **Meeus accuracy**.
Nutation does NOT enter — it is the mean↔true *equinox* difference (a frame this chain never
visits), not a positional error; it would matter only for a true-of-date / apparent-place
reduction (a separate, future chain carrying N(t)). The end-to-end test (§6.2) confirms this
empirically: the GCRS residual vs DE430 is ~3.7″ = the Meeus truncation, with **no** ~17″
nutation term. (An earlier draft of this note wrongly predicted a nutation-dominated residual;
the DE430 measurement corrected it — recorded here rather than quietly fixed.)

## 6. Independent oracle — element-wise (conformance) + DE430 (fidelity)

Per the INDEPENDENCE LESSON (erfa precession/obliquity ARE the IAU 2006 analytical theory, so
matching them only proves code-vs-theory), the gate `test_frame_chain` pins TWO levels:

  1. **Element-wise vs ERFA (conformance).** precession_iau2006(t) == `erfa.pmat06`,
     obliquity_rotation(ε) == the Rx(−ε) of `erfa.ecm06`, composite M(t) == `erfa.ecm06`ᵀ — to
     ~1e-12 per element over JD 2415020 / 2451545 / 2459000.5 / 2469807.5 (the §A epochs).
     This proves the rotations implement IAU 2006; it is NOT a fidelity claim.
  2. **End-to-end vs JPL DE430 (fidelity).** `moon_meeus` → (1) → M(t) → GCRS, compared to the
     Moon's GCRS position read DIRECTLY from the in-repo DE430 table (cols 10–12 ARE J2000
     mean-equatorial ≈ GCRS — the truly INDEPENDENT numerical ephemeris; NOT astropy GCRS
     frames, which blow up at t≠0). Because the chain is mean-consistent (§5) the orthonormal
     rotation adds ≤1e-12, so the angular residual is purely the Meeus truncation (Moon ≤ ~10″
     at 60 terms). Measured: **~3.7″**, and the gate bound (10″, the Meeus grade) MAJORIZES it
     (no perceived fidelity). body_position_gcrs is the consumer; the DE430 reference is
     born-digital from tools/gen_frame_oracle.py.

## 7. What falls out

  - `astronomy/frames.h` (extend): `obliquity_rotation(ε) = rot_x(−ε)`;
    `precession_iau2006(Epoch) → Matrix3` (FW polynomials via TrackedPolynomial + the §A fw2m
    composition of rot_x/rot_z); `ecliptic_to_gcrs(Epoch) → Matrix3 = precession_iau2006ᵀ ·
    obliquity_rotation(ε_A)`. Wires the `obliquity.h` orphan.
  - `ephemeris/body_position_gcrs.h` (new): compose a Meeus EclipticState (λ,β,Δ) through (1)
    and `ecliptic_to_gcrs` to a GCRS Cartesian position.
  - `tests/test_frame_chain` (new ExeGate): the §6 element-wise + DE430 oracle.

The L2-deferred precession/nutation/obliquity_rotation/Frame pieces (held back as anti-dead-
code) now have their consumer and land here.
