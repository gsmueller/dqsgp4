# Born-Digital Sources Catalog (research findings, 2026-06-07)

A curated, **verified** map of born-digital sources for the library's theory, per module. The trust model
([[feedback_born_digital_latex]], [[feedback_no_perceived_fidelity]]): coefficients come from a citable
born-digital source (published reference code / tables / LaTeX), and every set is **independently verified**
against a trusted oracle to a stated tolerance — never typed from memory, never asserted. Where this note says
"VERIFIED" the coefficients were re-derived locally and matched the oracle this turn (commands reproducible).

## Local oracles (independent, OFFLINE, already installed)

- **ERFA 2.0.1.5** (`pyerfa`, the IAU **SOFA** reference port, BSD) — exposes the COMPLETE frame chain and
  ephemeris helpers, confirmed present this turn: `obl06 pfw06 pmat06 nut06a nut00b num06a bi00 pnm06a gmst06
  gst06a ee06a era00 pom00 c2t06a ecm06 epv00 plan94 moon98 dtdb dat`. This is the verification oracle for
  L1/L2 and a cross-check for L3. Repos: `github.com/liberfa/erfa` (C source = the coefficients),
  `github.com/liberfa/pyerfa` (Python).
- **astropy 7.1.1** — `get_body` (JPL DE via `jplephem`) + coordinate frames; the L3 position oracle.
- **In-repo** (born-digital): Vallado `MSIS_Vers.cpp`/`EopSpw.cpp` (NRLMSISE-00 / EOP), JPL **DE430**
  `sgp4_references/vallado_celestrak/datalib/sunmooneph_430t12.txt`.

## L1 — Time scales & epochs

| component | born-digital source | local oracle |
|---|---|---|
| TAI−UTC leap seconds | SOFA/ERFA `dat.c`; IERS Bulletin C | `erfa.dat` |
| TDB−TT periodic | SOFA/ERFA `dtdb.c` (Fairhead & Bretagnon) | `erfa.dtdb` |
| UT1−UTC (ΔUT1) | IERS Bulletin A/B (finals2000A) | — (EOP file) |
| JD/epoch conventions | IERS Conventions 2010 Ch.1; astropy.time | erfa, astropy |

Primary: **IERS Conventions (2010)**, `iers-conventions.obspm.fr/content/chapter5/icc5.pdf` (and Ch.1–4).

## L2 — Reference frames — SOURCE: SOFA/ERFA (the IAU reference implementation)

| component | ERFA routine (= the coefficients) | primary literature | status |
|---|---|---|---|
| mean obliquity ε_A | `obl06.c` | IERS Conv. 2010 Eq.5.40; Capitaine+ 2003 | **VERIFIED** `obliquity.h` ≤1.1e-11″ vs `erfa.obl06` |
| precession (Fukushima-Williams γ̄,φ̄,ψ̄) | `pfw06.c`, matrix `pmat06.c`=`fw2m(pfw06)` | Capitaine, Wallace & Chapront (2003) A&A 412,567; IERS Conv.2010 §5.6 | **VERIFIED bit-exact** (0.00e+00″) vs `erfa.pfw06`/`pmat06` — coeffs §A below |
| nutation Δψ,Δε | `nut06a.c` (IAU2000A/2006), `nut00b.c` (2000B abridged, 77 terms) | IERS Conv.2010 §5.5; Mathews+ 2002 | `erfa.nut06a/nut00b` (local) |
| frame bias ICRS↔J2000 | `bi00.c` | IERS Conv.2010 §5.4 | `erfa.bi00` (local) |
| GMST / GAST / ERA | `gmst06.c`,`gst06a.c`,`era00.c` | IERS Conv.2010 §5.5; Aoki+ 1982 (our `sidereal_time.h`) | `erfa.gmst06` (local) |
| polar motion W(t) | `pom00.c` | IERS Conv.2010 §5.4 | `erfa.pom00` (local) |
| ecliptic-of-date → ICRS | `ecm06.c` | — | **VERIFIED** ecliptic→equatorial = `Rx(−ε)` exact at t=0 vs astropy |

The full GCRS→ITRS map `c2t06a` = `pom00·R3(GST)·pnm06a` is entirely ERFA, so every L2 increment is locally
gate-able element-wise against the IAU reference.

## L3 — Generic ephemeris

| body | low-precision (truncation) | full theory | local oracle | status |
|---|---|---|---|---|
| Sun | **Meeus §25** (pymeeus `Sun.py`; our `earth_orbit.h`) | VSOP87 (Bretagnon & Francou) | astropy `get_body`; `erfa.epv00` | **VERIFIED** `sun_meeus.h` ≤0.01° vs astropy |
| Moon | **Meeus §47** (pymeeus `Moon.py`, Tables 47.A/47.B) | ELP2000-82B (Chapront-Touzé & Chapront) | astropy `get_body`; **`erfa.moon98`** (in-ERFA Moon) | **VERIFIED** 12+8 terms ≤2.1′ vs astropy — coeffs §B below |
| planets | — | VSOP87 | `erfa.plan94` | — |

Born-digital homes:
- **Meeus, *Astronomical Algorithms* (2nd ed., 1998)** Ch.25 (Sun), Ch.47 (Moon) — realised VERBATIM by
  **pymeeus** `github.com/architest/pymeeus` (`pymeeus/Sun.py`, `pymeeus/Moon.py`). pymeeus is NOT installed,
  but its tables are public source; the §B coefficients below were fetched from it and locally verified.
- **VSOP87** (solar/planetary) — VizieR `VI/81`; **ELP2000-82B** (lunar) — VizieR `VI/79`
  (Chapront-Touzé & Chapront 1988); C port `github.com/variar/elp2000-82b`. NASA `eclipse.gsfc.nasa.gov`
  (`SEpath/ve82-predictions.html`) documents the VSOP87+ELP2000/82 realisation.
- **JPL DE430/DE440** — via astropy/`jplephem`; in-repo DE430 table (above).

## Appendix A — IAU 2006 precession FW angles (VERIFIED bit-exact vs erfa.pfw06/pmat06)

From ERFA `src/pfw06.c` (= SOFA `iauPfw06`); t = TT Julian centuries from J2000.0; result in arcsec, ×π/648000
→ rad. Reproduces `erfa.pfw06` AND `erfa.pmat06` to **0.00e+00″ / 0.00e+00** at JD 2415020/2451545/2459000.5/
2469807.5 (this turn). The precession-bias matrix is `fw2m(γ̄,φ̄,ψ̄,ε_A) = R1(−ε_A)·R3(−ψ̄)·R1(φ̄)·R3(γ̄)`.

```
γ̄(t) = ( -0.052928   + 10.556378·t   + 0.4932044·t²  - 0.00031238·t³ - 0.000002788·t⁴ + 0.0000000260·t⁵ )″
φ̄(t) = ( 84381.412819 - 46.811016·t   + 0.0511268·t²  + 0.00053289·t³ - 0.000000440·t⁴ - 0.0000000176·t⁵ )″
ψ̄(t) = ( -0.041775   + 5038.481484·t  + 1.5584175·t²  - 0.00018522·t³ - 0.000026452·t⁴ - 0.0000000148·t⁵ )″
ε_A(t) = obl06  (already in src/astronomy/obliquity.h, VERIFIED ERFA-exact)
```

Each c_k is a finite-digit IAU fit → `model_coefficient` (digits → accuracy, storage → T-precision), exactly as
`obliquity.h` already treats Eq.5.40. Implementation: `precession_iau2006(Epoch)→Matrix3` + the `fw2m`
composition (R1/R3 = our `rot_x`/`rot_z`), gated element-wise vs `erfa.pmat06`.

## Appendix B — Lunar Meeus §47 leading terms (VERIFIED ≤2.1′ vs astropy with 12+8 terms)

From pymeeus `Moon.py` (Meeus Tables 47.A longitude/distance, 47.B latitude). Mean elements (deg, T in Julian
centuries from J2000): L′,D,M,M′,F polynomials per Meeus §47 (fetched verbatim — recorded in the verification
script). Each periodic row = integer multipliers (D,M,M′,F) and a coefficient; terms with the Sun's anomaly M
carry the eccentricity factor E^|M|, E = 1 − 0.002516·T − 0.0000074·T². λ = L′ + Σl·1e-6°, β = Σb·1e-6°,
Δ = 385000.56 km + Σr·1e-3 km.

```
Table 47.A (D, M, M′, F | Σl[1e-6°] | Σr[1e-3 km]) — top 12:
 0 0 1 0 | 6288774 | -20905355     (equation of centre, main)
 2 0 -1 0| 1274027 |  -3699111     (EVECTION)
 2 0 0 0 |  658314 |  -2955968     (VARIATION)
 0 0 2 0 |  213618 |   -569925
 0 1 0 0 | -185116 |     48888     (ANNUAL EQUATION, ×E)
 0 0 0 2 | -114332 |     -3149
 2 0 -2 0|   58793 |    246158
 2 -1 -1 0|  57066 |   -152138     (×E)
 2 0 1 0 |   53322 |   -170733
 2 -1 0 0|   45758 |   -204586     (×E)
 0 1 -1 0|  -40923 |   -129620     (×E)
 1 0 0 0 |  -34720 |    108743
Table 47.B (D, M, M′, F | Σb[1e-6°]) — top 8:
 0 0 0 1 | 5128122   (main, ≈ i·sin F)
 0 0 1 1 |  280602
 0 0 1 -1|  277693
 2 0 0 -1|  173237
 2 0 -1 1|   55413     (×E)
 2 0 -1 -1|  46271     (×E)
 2 0 0 1 |   32573
 0 0 2 1 |   17198
```

Full tables (60+60 terms) in pymeeus `Moon.py` → reach ~10″; the truncation accuracy is Σ omitted |amplitude|
(the generator pattern), gate-able to shrink as terms are added, vs astropy / `erfa.moon98`.

## Reproduce

`tools/gen_solar_oracle.py` (solar); the precession + lunar verification this turn is in
[[session-2026-06-06b-library-rearchitecture]]. Re-run vs `erfa`/`astropy` to confirm before each implementation.
