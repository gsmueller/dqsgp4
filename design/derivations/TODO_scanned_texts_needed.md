# TODO: Scanned Texts Requiring Markdown Conversion

These source documents are currently available only as scanned image PDFs (from ADS,
library archives, or similar). The next development phase requires working from the
actual mathematical content. OCR from scanned PDFs is error-prone; for critical
derivations, a clean markdown transcription provides a verified, searchable, citable
reference.

**Priority classification:**
- **CRITICAL** — coefficients in the textbook are derived directly from this source;
  every equation must be accessible and verifiable
- **HIGH** — significant derivations reference this source; errors would affect multiple
  chapters
- **MEDIUM** — background reference or cross-check; primary derivation uses a clean
  digital source

---

## 1. Brouwer (1959) — CRITICAL

**Full citation:** Brouwer, D. (1959). "Solution of the Problem of Artificial Satellite
Theory Without Drag." *Astronomical Journal*, 64, 378–397.

**ADS location:** `https://ui.adsabs.harvard.edu/abs/1959AJ.....64..378B`

**Chapters using this source:** Ch 16 (J₂ secular first order), Ch 17 (second-order
secular), Ch 18 (short-period corrections), Ch 19 (long-period corrections), Ch 20
(osculating elements), Ch 33 (secular update)

**Why markdown is needed:** The secular rate coefficients (k₂, k₄, a₁, d₁, d₂, etc.)
and the short-period correction polynomials are all transcribed from this paper. The
scanned ADS PDF has OCR errors in the coefficient tables. Lara (2021) corrects several
of these but does not reproduce the full set of Brouwer's original intermediate
quantities. A verified markdown transcription of §§1–5 of Brouwer (1959) is needed for
the Develop phase.

**Sections needed:** §1 (Hamiltonian setup), §2 (secular terms), §3 (short-period
terms), §4 (long-period terms), §5 (final element formulas). The coefficient tables in
§§2–5 are the critical content.

**Known issues:** Lara (2021, J. Astronaut. Sci. 68) identifies sign errors in
Brouwer's long-period inclination rate formula. Any markdown transcription must flag
these and document the correction.

---

## 2. Hoots & Roehrich (1980) / SR3 — CRITICAL

**Full citation:** Hoots, F. R. & Roehrich, R. L. (1980). "Models for Propagation of
NORAD Element Sets." *Spacetrack Report No. 3* (SR3). Aerospace Defense Command,
Peterson AFB, Colorado.

**Online location:** `https://celestrak.org/SPACETRACK/P3Doc/` (widely hosted;
digitized from the original DOD technical report)

**Chapters using this source:** Ch 31 (TLE parsing), Ch 32 (element recovery), Ch 33
(secular update), Ch 34 (near-space pipeline), Ch 35 (deep-space pipeline), Ch 37
(precomputed constants), App A (WGS72 constants)

**Why markdown is needed:** SR3 is the primary SGP4 implementation reference. However,
it contains known transcription errors documented by Vallado et al. (2006). A markdown
version allows side-by-side comparison with the Vallado corrections. Many SR3 formulas
are used in the textbook as *falsifiable hypotheses* that must be independently derived;
the markdown version is the primary record of what SR3 actually claims.

**Sections needed:** Full document — approximately 90 pages. Priority: §§4–7 (the SGP4
model equations), Appendix A (deep space), Appendix B (DPINIT/DPSEC/DPPER equations).

**Known issues:** Vallado, Crawford, Hujsak & Kelso (2006), "Revisiting Spacetrack
Report #3," AIAA 2006-6753, corrects multiple formulas. The markdown must annotate each
corrected formula with the Vallado correction and its AIAA paper location.

---

## 3. Lane & Hoots (1979) — HIGH

**Full citation:** Lane, M. H. & Hoots, F. R. (1979). "General Perturbation Theories
Derived from the 1965 Lane Drag Theory." *Aerospace Defense Command Technical Report.*
Project Space Track Report.

**Location:** May be available through DTIC or NRL technical report archives. Check
DTIC accession number ADA081960 (tentative).

**Chapters using this source:** Ch 21 (atmospheric density model), Ch 22 (drag
coefficients), specifically the Lane combinatorial method for the orbit-averaged drag
integral.

**Why markdown is needed:** The binomial series × Wallis integral formulation for the
C₁–C₅ drag coefficients (the Lane combinatorial method) is derived in this report.
Lane's 1965 original has known errors corrected by Lane & Hoots (1979). A markdown
transcription of the corrected derivation is the authoritative source for Ch 22.

**Sections needed:** The section deriving the drag integral expansion (typically §2–3
of such reports), the C₁–C₅ coefficient formulas, and the density model power-law
assumption.

**Known issues:** The 1965 Lane paper has errors in the coefficient formulas; the 1979
Lane & Hoots paper provides corrections. The markdown must be of the 1979 (corrected)
version, not the 1965 original.

---

## 4. Kaula (1966) — HIGH

**Full citation:** Kaula, W. M. (1966). *Theory of Satellite Geodesy*. Blaisdell
Publishing Company, Waltham, Massachusetts. (Reprinted by Dover, 2000.)

**Chapters using this source:** Ch 13 (geopotential expansion), Ch 15 (Kaula expansion
and F/G functions), Ch 16 (leading into Ch 16 from Ch 15 Kaula functions)

**Why markdown is needed:** The Fₗₘₚ (inclination functions) and Gₗₚq (eccentricity
functions) tables in Kaula (1966) §3.3 are the foundation for Ch 15. The Dover reprint
is a photographic reproduction of the 1966 typeset edition; many formulas in §3 are
partially illegible in low-resolution scans. A markdown transcription of Ch 3 (§§3.1–
3.4) is needed.

**Sections needed:** Chapter 3 (§§3.1–3.4: the geopotential disturbing function, Kaula
F and G functions, the resonance condition). Table 3.1 (Fₗₘₚ for l=2,3,4) and Table
3.2 (Gₗₚq for e up to e⁶) are the critical content.

**Known issues:** The Sneeuw (2022) geodesy lecture notes (clean digital source,
already in the project) provide a modern re-derivation of the Kaula functions that
should be used as the primary cross-check. Kaula (1966) should be used to verify
the original definitions and notation.

---

## 5. Heiskanen & Moritz (1967) — MEDIUM

**Full citation:** Heiskanen, W. A. & Moritz, H. (1967). *Physical Geodesy*. W. H.
Freeman and Company, San Francisco.

**Chapters using this source:** Ch 14 (equipotential ellipsoid, Somigliana formula,
gravity potential), indirectly Ch 13

**Why markdown is needed:** Ch 14 derives the equipotential ellipsoid and Somigliana
formula from first principles. The Heiskanen & Moritz derivation (§2.7–2.9) is the
canonical reference for this. The book is available only as a photographic scan.

**Sections needed:** §§2.7–2.9 (the normal gravity formula and Somigliana derivation),
§3.1 (spherical harmonic expansion of the normal potential).

**Known issues:** The Sneeuw (2022) notes provide a modern version of this material in
clean digital form. Heiskanen & Moritz is the original, and spot-checking the Sneeuw
notes against it is part of the verification procedure. Only §§2.7–2.9 are strictly
needed; the rest of the book is background.

---

## 6. Aoki et al. (1982) — MEDIUM

**Full citation:** Aoki, S., Guinot, B., Kaplan, G. H., Kinoshita, H., McCarthy, D. D.,
& Seidelmann, P. K. (1982). "The new definition of Universal Time." *Astronomy &
Astrophysics*, 105, 359–361.

**ADS location:** `https://ui.adsabs.harvard.edu/abs/1982A%26A...105..359A`

**Chapters using this source:** Ch 29 (sidereal time, GMST polynomial)

**Why markdown is needed:** The GMST polynomial coefficients (Definition 29.3.1 in the
textbook plan) come directly from Eq. (14) of this paper. The ADS version may be
available as a scanned journal page. The four polynomial coefficients must be
transcribed exactly — any OCR or manual transcription error directly corrupts the
sidereal time computation.

**Sections needed:** The abstract, §3 (the new GMST formula), and Eq. (14) specifically.
This is a short paper (3 pages); the entire text should be transcribed.

**Known issues:** The IAU 2006 Earth rotation model (Capitaine et al. 2003) supersedes
this formula for precise work but it remains the SGP4 matched-pair formula. See
[A.29.2] in Ch 29 plan. The textbook must use the 1982 formula for SGP4 compatibility
and note the IAU 2006 alternative.

---

## 7. Vallado et al. (2006) — HIGH (corrections document, not primary source)

**Full citation:** Vallado, D. A., Crawford, P., Hujsak, R., & Kelso, T. S. (2006).
"Revisiting Spacetrack Report #3." AIAA 2006-6753, presented at AIAA/AAS
Astrodynamics Specialist Conference.

**Location:** `https://celestrak.org/publications/AIAA/2006-6753/` (freely available
digital PDF from Celestrak — this is a CLEAN DIGITAL SOURCE, no markdown needed)

**Note:** This document is already available as a clean digital PDF. It does NOT need
markdown conversion. However, it must be read alongside SR3 (item 2 above) to identify
all corrections. For each corrected formula in SR3, the textbook must show the original
SR3 formula, the Vallado correction, and the independently derived textbook result.

---

## Summary Table

| Priority | Document | Approx. Pages to Transcribe | Critical Content |
|----------|----------|-----------------------------|-----------------|
| CRITICAL | Brouwer (1959) | ~20 pp (§§1–5) | Secular rate coefficients, correction tables |
| CRITICAL | Hoots & Roehrich (1980) SR3 | ~90 pp (full) | SGP4 model equations, DPINIT/DPSEC |
| HIGH | Lane & Hoots (1979) | ~30 pp (§§2–3) | C₁–C₅ drag coefficient formulas |
| HIGH | Kaula (1966) Ch 3 | ~40 pp (§§3.1–3.4) | Fₗₘₚ, Gₗₚq tables and formulas |
| MEDIUM | Heiskanen & Moritz (1967) §§2.7–2.9 | ~15 pp | Somigliana formula derivation |
| MEDIUM | Aoki et al. (1982) | ~3 pp (full) | GMST polynomial Eq. (14) |

**Total estimate:** ~200 pages of critical scanned content requiring verified markdown
transcription before the Develop phase can be completed for Chapters 15–22 and 29–35.
