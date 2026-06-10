# Draft Plan: Chapter 31 — TLE Format and Parsing

**Part IX: The SGP4 Propagator** | Implementation file: `tle_parser.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $n_o$ | TLE mean motion (rev/day) | §31.2 |
| $e_o$ | TLE eccentricity (decimal point implied) | §31.2 |
| $i_o$ | TLE inclination (degrees) | §31.2 |
| $\Omega_o$ | TLE right ascension of ascending node (degrees) | §31.2 |
| $\omega_o$ | TLE argument of perigee (degrees) | §31.2 |
| $M_o$ | TLE mean anomaly (degrees) | §31.2 |
| $B^*$ | TLE drag-like coefficient (Hoots and Roehrich 1980 notation) | §31.2 |
| $\mathcal{V}$ | Tracked value $(v, \sigma_m, \delta_p, \delta_a)$ (Ch 1, Def. 1.2.4) | Ch 1 |
| $d_{\mathrm{checksum}}$ | Mod-10 checksum digit | §31.3 |

---

## Objectives

1. Specify the Two-Line Element format completely: all fields, column positions, units, and encoding conventions.
2. Derive the mod-10 checksum algorithm and prove its error-detection properties.
3. Convert each parsed field to a TrackedValue with correct $\sigma_m$ (measurement error from digit count).
4. Classify input errors: format errors, checksum failures, out-of-range values, epoch plausibility checks.
5. Establish the Tier classification of TLE-derived elements.

## Section Structure

### §31.1 Introduction

This section establishes the TLE as the primary input to the SGP4 propagator and provides a road map of the chapter.

Stub: The TLE is the primary input to the SGP4 propagator. It encodes mean orbital elements under the SGP4 model conventions. The format is fixed-column ASCII; parsing must be exact. Forward-reference to Ch 3 (matched pair: TLE elements are meaningful only within SGP4) and Ch 32 (element recovery from the parsed values).

### §31.2 TLE Field Specification

This section provides a complete field-by-field specification of the TLE format and assigns measurement error bounds to each field.

Stub: Table 31.2.1: complete field-by-field specification of Line 1 and Line 2. Columns, name, unit, encoding (decimal, integer, assumed decimal point, exponent notation for $B^*$). All field widths given as exact column ranges per the Space-Track standard. [M.31.1] Inclination resolved to 0.0001°; $\sigma_m \approx 5\times10^{-5}$°. [M.31.2] Mean motion resolved to $10^{-8}$ rev/day. [M.31.3] $B^*$ encoded as mantissa × $10^{\text{exp}}$; precision limited to 5 significant digits. Definition 31.2.1 (assumed decimal point): eccentricity field contains 7 digits with an implied leading "0."; value = integer ÷ $10^7$.

### §31.3 Checksum Algorithm

This section defines the mod-10 checksum, proves its single-error detection properties, and states its limitations.

Stub: Definition 31.3.1 (mod-10 checksum): sum of all integer digits in the line (letters contribute 0, minus signs contribute 1, spaces/periods contribute 0); result mod 10 equals the check digit. Proposition 31.3.1: detects all single-digit substitution errors and all single transposition errors of adjacent non-equal digits. — *Proof approach: for single-digit substitution, the checksum changes by $(d' - d) \bmod 10 \ne 0$ for $d' \ne d$; for adjacent transposition of digits $a \ne b$, the checksum changes by $(b - a) + (a - b) = 0$ only if $a = b$, so all $a \ne b$ transpositions are detected.* [A.31.1] Does not detect all two-digit errors; a TLE with two compensating errors can pass the checksum.

### §31.4 Conversion to TrackedValues

This section converts each parsed TLE field to SI units and wraps it in a TrackedValue with measurement error $\sigma_m$ derived from the field digit count.

Stub: Each parsed field is converted to SI units and wrapped in a TrackedValue. Procedure: (1) parse integer or decimal string to double; (2) assign $\sigma_m$ = half the least significant digit in the original unit; (3) convert unit (degrees to radians using Ch 7, exact rational factor); (4) assign $\delta_p$ = unit conversion rounding; (5) $\delta_a = 0$ at parse time (model error assigned at Ch 32). Table 31.4.1: $(v, \sigma_m, \delta_p, \delta_a)$ for each element at parse time.

### §31.5 Input Error Classification

This section classifies parsing failures into four error classes and specifies the propagator's response to each.

**Example 31.5.1** (ISS TLE parse): Parse the ISS TLE (NORAD 25544, example epoch 2023-001.50000000): verify checksum on both lines, extract all fields, convert $n_o = 15.49196328$ rev/day to $1.13170\times10^{-3}$ rad/s, assign $\sigma_m = 5\times10^{-9}$ rev/day from digit count. Source: Celestrak SGP4 verification TLE set (Vallado et al. 2006).

Stub: Three error classes. Class F (format error): column layout violation, non-numeric character in numeric field, missing line, line length error. Class C (checksum error): computed checksum ≠ stored digit; suspect transcription. Class R (range error): $e_o \notin [0, 1)$, $i_o \notin [0°, 180°]$, $n_o \leq 0$, epoch before 1957. Class P (plausibility warning): $n_o > 17$ rev/day (below ~150 km orbit would decay rapidly), $B^* < 0$ (non-physical). Each class triggers a distinct error code; the propagator refuses to proceed on Class F and C; Class R and P are caller-configurable.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1, Theorem 1.2.4 | §31.4 | TrackedValue definition for wrapping parsed fields |
| Ch 3, Theorem 3.2.1 | §31.1 | Matched pair principle: TLE elements tied to SGP4 model |
| Ch 7, angle arithmetic | §31.4 | Degree-to-radian conversion with exact rational factor |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 32, element recovery | §31.2, §31.4 | Parsed TLE fields ($n_o, e_o, i_o, B^*$) with TrackedValues |
| Appendix A | §31.2 | WGS72 constants that define TLE unit system |
| Appendix C | §31.2 | tle_parser.h code mapping |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.31.1] | M | §31.2 | Inclination precision from TLE digit count |
| [M.31.2] | M | §31.2 | Mean motion precision from TLE digit count |
| [M.31.3] | M | §31.2 | $B^*$ limited to 5 significant digits |
| [A.31.1] | A | §31.3 | Mod-10 checksum does not detect all two-digit errors |
| [A.31.2] | A | §31.2 | TLE mean motion $n_o$ is in the Kozai (1959) convention, not the Brouwer (1959) convention; the two frameworks define mean elements differently and give numerically distinct mean motions; Ch 32 must perform the Kozai→Brouwer conversion before applying SGP4 secular rates; mixing conventions introduces a systematic semi-major axis error of $O(J_2)$ |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 4 |
| Theorems | 0 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 2 |
| Examples | 1 |
| Error Notes | 5 |
| Equations | ~5 |
| Sections | 5 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §31.1 | Draft | Introduction |
| §31.2 | Draft | TLE field specification |
| §31.3 | Draft | Checksum algorithm |
| §31.4 | Draft | Conversion to TrackedValues |
| §31.5 | Draft | Input error classification |
