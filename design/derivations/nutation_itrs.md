# Nutation (IAU 2000A luni-solar) and the Earth-fixed chain GCRS → ITRS

**Theory note (FIRST, before code) for replan item R4b** (`design/PROFESSIONAL_LIBRARY_PLAN.md` §6) — the
Earth-fixed output chain: nutation N(t), the IAU-2006-era sidereal angles (ERA/GMST/GAST), polar motion W(t),
and the full GCRS→ITRS rotation. This is where nutation lands **with its honest consumer** (Earth-fixed
outputs/ground tracks) — the anti-dead-code resolution recorded at the replan. Governed by
[[feedback_theory_first_library]] + [[feedback_no_perceived_fidelity]]. Companion code:
`astronomy/nutation_terms_iau2000a.h` (auto-generated), `astronomy/nutation.h`, `astronomy/earth_rotation.h`;
oracle `tools/gen_itrs_oracle.py` (erfa nut06a/gmst06/gst06a/pom00/c2t06a + the in-repo CSSI EOP file); gate
`tests/test_nutation` (NUT1).

---

## 1. The physics — why the Earth's spin axis nods

The Sun and Moon exert a torque on the Earth's equatorial bulge (the same J₂ oblateness the geopotential
carries). The secular part of the response is **precession** (the equinox sweeping ~50″/yr — already modeled,
`precession_iau2006`, erfa-gated). The **periodic** part is **nutation**: the pole nods as the lunisolar
torque geometry cycles through the orbital configurations. Because the torque depends on the Sun/Moon
positions, the response is a **Poisson series in the five Delaunay arguments**

```
  l   Moon's mean anomaly       l′  Sun's mean anomaly      F = L_moon − Ω
  D   mean elongation           Ω   Moon's node             (the same arguments as moon_meeus §47)
```

expressed as the nutation in longitude Δψ (along the ecliptic) and in obliquity Δε:

```
  Δψ = Σ_k (A_k + A′_k·t)·sin θ_k + A″_k·cos θ_k        θ_k = n₁l + n₂l′ + n₃F + n₄D + n₅Ω
  Δε = Σ_k (B_k + B′_k·t)·cos θ_k + B″_k·sin θ_k
```

The dominant term is the 18.6-year nodal nutation: −17.2064161″ sin Ω in Δψ, +9.2052331″ cos Ω in Δε.

## 2. The model — IAU 2000A luni-solar, generatively truncatable

- **Source (born-digital, in-repo):** the SOFA `nut00a` luni-solar table,
  `sgp4_references/vallado_celestrak/datalib/iau00ansofa.dat` — **678 terms**, coefficients **exact integers
  in 1e-7 arcsec** (so amplitude accuracy is pure truncation, the lunar-table pattern). Transcribed
  mechanically by `tools/gen_nutation_terms.py`, **amplitude-sorted** so `n_terms` dials fidelity up with the
  tracked tail Σ_{k>n}|A_k| — the third instance of the truncated-series-with-tracked-tail pattern
  (kepler_series, poisson_series, now nutation).
- **The planetary series is BOUNDED, not modeled:** SOFA's 687 planetary terms
  (`iau00anplsofa.dat`) total |amplitude| **Δψ ≤ 55379e-7″ = 5.54 mas, Δε ≤ 19022e-7″ = 1.90 mas** (measured
  by the generator; the conservative triangle bound). These floors are deposited in the angles' accuracy
  channel. (The actual planetary signal is a few tenths of a mas — the gate measures it; the floor majorizes.)
- **Fundamental arguments:** the five IERS quartic polynomials [arcsec], transcribed **verbatim from the
  in-repo AstroLib `fundarg` (e06 branch)** (l: 485868.249036 + 1717915923.2178t + …; the canonical
  Simon-1994/IERS-2010 values), evaluated as `TrackedPolynomial`s.
- **The matrix:** N(t) = R1(−(ε_A+Δε)) · R3(−Δψ) · R1(ε_A), with ε_A from `obliquity_iau2006` (SC1). Our
  rot_x/rot_z are exactly SOFA's R1/R3 (proved bit-exact in the FRAME2 precession gate), so the composition
  is the erfa `numat` form.

## 3. The rotation angles — ERA, GMST06, GAST06

The IAU-2006-era chain (NOT the in-repo Aoki-82 GMST, which serves the frozen SGP4 TEME convention):

- **ERA** (Earth rotation angle): the IAU 2000 Res. B1.8 **defining relation**
  `ERA = 2π(0.7790572732640 + 1.00273781191135448·(JD_UT1 − 2451545))` — both constants exact-by-convention
  (`defined`, CR1B-allowlisted). UT1 is the caller's (ΔUT1 from EOP).
- **GMST06** = ERA + the polynomial (0.014506 + 4612.156534t + 1.3915817t² − 0.00000044t³ + 0.000029956t⁴ +
  0.0000000368t⁵)″, transcribed verbatim from the in-repo AstroLib `gstime00`. **Transcription arbitrated by
  measurement:** this form vs `erfa.gmst06` differs by ≤ **2.4 μas** over 1995–2024 (the t⁴/t⁵ sign question
  is immaterial at the chain's mas grade; recorded, bounded).
- **GAST06** = GMST06 + the equation of the equinoxes: Δψ·cos ε_A + the two classical complementary terms
  **+0.00264″ sin Ω + 0.000063″ sin 2Ω** (in-repo AstroLib verbatim). The remaining ~31 complementary
  (`eect00`) terms total ~30 μas — omitted with a 50 μas accuracy floor.

## 4. Polar motion and the full chain

- **W(t) = R3(−s′)·R2(x_p)·R1(y_p)** (IERS 2010 §5.3) with the TIO locator s′ = −47 μas·t
  (`model_coefficient`). x_p, y_p are **measured IERS data, not a model** — caller-supplied (the erfa
  `pom00` signature pattern, same reasoning as the MSIS space-weather decision D2-a); the in-repo CSSI EOP
  file (`EOP-All-v1.1_2025-01-10.txt`) feeds the *gate* its born-digital values.
- **The full chain:** `r_ITRS = W(t) · R3(GAST) · N(t) · P(t) · r_GCRS` — composed of the four gated pieces
  (P from FRAME2; N, GAST, W here). Equinox-based; equals the CIO-based erfa `c2t06a` to the model's own
  grade (the IERS Ch. 5 equivalence), which is exactly what the gate measures.

## 5. Error budget (the honest grade)

| piece | bound (deposited in accuracy) |
|---|---|
| Δψ truncation (n_terms < 678) | Σ_{k>n} (|A|+|A′t|+|A″|) · 1e-7″ — the dial-up tail |
| Δψ omitted planetary | 5.54 mas floor (measured Σ|A|) |
| Δε ditto | tail + 1.90 mas floor |
| GMST06 transcription | 2.4 μas (measured vs erfa) |
| Eqeq omitted complementary | 50 μas floor (~30 μas measured total) |
| s′, EOP | s′ digit floor; x_p/y_p carry the caller's (IERS) uncertainty |

Chain grade: **~mas-level (≈ 5e-8 rad in the matrix elements), dominated by the honest planetary floor** —
NOT machine precision, and stated so. (The *full* 2000A incl. planetary terms, if a consumer ever needs
sub-mas, is a data-file extension of the same generator/engine — the dial-up pattern continues.)

## 6. The gate (NUT1) — erfa arbitration at every layer

1. **Angles:** Δψ, Δε (full 678 terms) vs `erfa.nut06a` (the COMPLETE 2000A+planetary+P03 model) over the 5
   oracle epochs — measured residual = the real planetary signal (~sub-mas), the tracked floors MAJORIZE it.
2. **Dial-up:** residual(n=20) ≫ residual(678); the tracked truncation tail majorizes the n=20 deficit.
3. **GAST** vs `erfa.gst06a` — sub-mas (the Δψ planetary part × cos ε + the 50 μas eect floor majorize).
4. **W** vs `erfa.pom00` — element-wise at machine precision (pure exact-angle rotations).
5. **The full chain** vs `erfa.c2t06a` — element-wise ≤ the chain grade; the composed tracked accuracy
   majorizes the measured residual. (Real EOP from the in-repo file, identical values to the oracle's.)
6. Matrix identities (N, W, chain orthonormal; det +1) + precision tightens with wider T.

The retro-benefit recorded at the replan: the third-body/SRP TEME frame bound (~30″ for omitted
nutation+Eqeq) can now be *computed* rather than bounded — left as a future tightening; the bound is honest.
