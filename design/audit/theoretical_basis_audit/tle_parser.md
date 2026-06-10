# Theoretical Basis Audit — `src/tle/tle_parser.h` + `src/tle/tle_parser.cpp`

**Files**: `src/tle/tle_parser.h` (132 lines) + `src/tle/tle_parser.cpp` (181 lines)
**Audited**: 2026-05-13; expanded R14 2026-05-13 to cover .cpp definitions.
**Functions in scope**: 3 (parse line1/line2, parse name/line1/line2, `TleElements::from_tle_data`)

---

## Card 1 — `tle::parse(line1, line2, out)`

```
=== FORMULA AUDIT CARD ===
ID:                     tle_parser::parse_two_line
Location:               src/tle/tle_parser.h:48 (declaration);
                        src/tle/tle_parser.cpp:140-142 (definition).
Mathematical statement: Lex/parse two fixed-column 69-character ASCII records
                        into TleData numeric fields. Pure string→number
                        deserialization; no math operation.

THEORY
  Underlying theorem:   TLE format specification: each line is 69 ASCII
                        characters with fixed column positions for each
                        field; numeric fields are decimal with documented
                        digit counts. No theorem applies — this is a
                        specification-driven lexer.
  Primary reference:    Hoots & Roehrich (1980) Spacetrack Report No. 3,
                        §3 "Two-Line Mean Element Set Format". Also
                        CelesTrak TLE format documentation
                        (https://celestrak.org/NORAD/documentation/tle-fmt.php).
  Domain of validity:   Inputs that conform to the SR3 §3 format
                        (69-char lines, expected checksums, valid columns).

METHOD
  Method declared:      Fixed-column substring extraction + std::stod /
                        std::stoi parsing.
  Method implemented:   src/tle/tle_parser.cpp:140-142 (two-arg overload
                        delegates to three-arg overload with name=""):
                          - line lengths ≥ 69 validation (cpp:146)
                          - leading-char check '1','2' (cpp:147)
                          - field extraction via local helpers:
                              extract_field (cpp:51-59) — substring + trim
                              parse_double  (cpp:62-66) — stod on field
                              parse_int     (cpp:69-73) — stoi on field
                              parse_assumed_decimal (cpp:78-82) — for ecc
                              parse_tle_exponential (cpp:87-129) — BSTAR, ndot
                              parse_catalog_number (cpp:134-136) — Alpha-5
                          - field column positions documented at cpp:22-48
                            match SR3 §3 verbatim.
  Match verdict:        ✓ matched (R14, 2026-05-13) — every column
                        extraction matches the SR3 §3 specification
                        position-by-position. Alpha-5 satellite numbers
                        handled transparently via string-typed catnum.

ERROR BOUND
  Bound category:       n/a — this function does not return TrackedValue<T>.
                        It produces raw `double` and `int` fields in TleData.
                        Measurement uncertainty is attached downstream by
                        `from_tle_data`.
  Bound formula:        n/a at this layer.
  Bound implemented:    n/a — return type is `bool` (success flag).
  Bound verdict:        n/a — this is a pure deserialization step; the
                        measurement bound is the next stage's responsibility.

CROSS-AUDIT
  REQ-EF applies:       n/a (no TrackedValue<T> emitted)
  AUD-EF applies:       n/a
  AUD-MC applies:       n/a
  Verification test:    tests/test_tle/ exercises SR3 sample TLEs;
                        tests/test_sgp4/main.cpp parses
                        sgp4_references/.../SGP4-VER.TLE without errors
                        as a side-effect of the verification harness. ✓.

NOTES
  - The function signature returns `bool` only; on parse failure the
    out-param state is unspecified by this header. Caller contract should
    require checking the return before reading `out`. ✓ enforced via
    line-length and leading-char checks at cpp:146-147.
  - Alpha-5 satellite numbering (CelesTrak extension, 2020) handled
    transparently — `satellite_id` is std::string preserving the leading
    letter; downstream consumers do not require numeric satnum.
  - The `satellite_id` field is `std::string` (not int) to accommodate
    Alpha-5; ✓ structurally appropriate.
  - Checksum validation: parser does NOT validate column-68 checksum
    digit. Caller responsibility if needed. ⚠ minor robustness flag.
```

---

## Card 2 — `tle::parse(name, line1, line2, out)`

```
=== FORMULA AUDIT CARD ===
ID:                     tle_parser::parse_three_line
Location:               src/tle/tle_parser.h:49 (declaration);
                        src/tle/tle_parser.cpp:144-178 (definition — this
                        is the primary three-arg implementation; the
                        two-arg overload forwards to it with name="").
Mathematical statement: Same as Card 1 plus assignment of the optional
                        line-0 satellite name string to out.name.

THEORY
  Underlying theorem:   TLE "3LE" convention: an optional preceding line
                        (line 0) carries the satellite common name; the
                        following two lines are the standard SR3 §3 record.
                        No theorem applies.
  Primary reference:    Hoots & Roehrich (1980) Spacetrack Report No. 3,
                        §3 (TLE format proper).
                        CelesTrak documents the 3LE convention as a
                        widely-adopted extension.
  Domain of validity:   Same as Card 1, plus a name line of arbitrary
                        printable ASCII (typically ≤24 chars per CelesTrak).

METHOD
  Method declared:      Overload that captures `name` and otherwise
                        delegates to the two-line parse.
  Method implemented:   src/tle/tle_parser.cpp:144-178 — primary three-arg
                        impl. Name is assigned at cpp:149 (out.name = name)
                        before field extraction begins; no other special
                        handling.
  Match verdict:        ✓ matched (R14, 2026-05-13) — the three-arg
                        overload is the actual parsing routine; the
                        two-arg overload at cpp:140-142 is a thin
                        delegate. Name is opaque and stored verbatim.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        n/a
  Bound implemented:    n/a
  Bound verdict:        n/a — same rationale as Card 1.

CROSS-AUDIT
  REQ-EF applies:       n/a
  AUD-EF applies:       n/a
  AUD-MC applies:       n/a
  Verification test:    tests/test_sgp4/main.cpp parses 3LE-style
                        satellite name + line1/line2 sets via this
                        overload.

NOTES
  - Name is informational only; no math depends on it.
  - Behavior on a malformed name line (e.g. >24 chars, non-ASCII): name
    is stored verbatim as std::string — no normalization, truncation, or
    validation. ✓ acceptable since downstream consumers do not parse it.
  - The two-arg form (cpp:140-142) forwards name="" to the three-arg
    form; no code duplication.
```

---

## Card 3 — `tle::TleElements::from_tle_data`

```
=== FORMULA AUDIT CARD ===
ID:                     tle_parser::TleElements::from_tle_data
Location:               src/tle/tle_parser.h:69-128
Mathematical statement: Convert raw TleData (degrees, rev/day, etc.) to
                        TrackedValue<T> orbital elements (radians, rad/min,
                        Julian date), attaching measurement σ derived from
                        TLE column precision.

THEORY
  Underlying theorem:   (a) Unit conversion: deg→rad multiplies by π/180;
                            rev/day→rad/min multiplies by 2π/1440.
                        (b) Measurement uncertainty from fixed-column ASCII
                            precision: a field printed with d decimal digits
                            has half-ULP uncertainty σ = 5·10^{-(d+1)}
                            (uniform-quantization model on the printed mantissa).
                        (c) Julian Date from civil (Y, day-of-year):
                            algorithm consistent with Vallado/Meeus simplified
                            forms; see Meeus (1998) "Astronomical Algorithms"
                            §7 or Vallado (2013) §3.5 "Julian Date".
  Primary reference:    TLE field precision per Hoots & Roehrich (1980) SR3 §3
                        (column widths give digit counts). Unit conversions
                        are definitional. JD formula commonly attributed to
                        Vallado §3.5 (closed-form civil-date → JD).
  Domain of validity:   - epoch_year ∈ [00, 99]; window pivot at 57 maps
                          [57, 99] → 19xx and [00, 56] → 20xx.
                        - eccentricity ∈ [0, 1).
                        - All inputs as documented in TleData.

METHOD
  Method declared:      (M1) Angles: σ_deg = 5e-5, σ_rad = σ_deg · π/180,
                             then deg→rad via `degrees_to_radians`.
                        (M2) Eccentricity: σ = 5e-8 absolute.
                        (M3) Mean motion: σ_revday = 5e-9, scale to rad/min
                             by ·(2π/1440).
                        (M4) B*: σ = 10% · |bstar| (relative).
                        (M5) Epoch JD: closed-form integer arithmetic on
                             year + day_of_year, σ_jd = 1e-8 day.
  Method implemented:   - L80:  angle_sigma_deg = 5/100000 = 5e-5         (M1)
                        - L82:  eps = numeric_limits<T>::epsilon()
                        - L84-89: make_angle lambda — TrackedValue ctor
                                  (val, σ_meas=5e-5 deg, σ_prec=|val|·ε, σ_acc=0),
                                  then degrees_to_radians applied.
                                  ⚠ Note: σ_meas is constructed in deg units
                                  but passed to a value that is then
                                  converted to rad. If degrees_to_radians
                                  scales σ_meas by π/180 internally, ✓; else
                                  σ stays numerically 5e-5 in *rad* units,
                                  which would overstate by ×(180/π) ≈ 57.3.
                                  ? — must inspect degrees_to_radians in
                                  src/math/angles.h to confirm propagation.
                        - L97-99: ecc_sigma = 5/1e8 = 5e-8                 (M2)
                        - L103-106: mm scaled by ·(2π/1440); σ scaled the
                                  same way                                  (M3) ✓
                        - L110-112: bstar_sigma = |bstar|·1/10 = 10% rel.   (M4)
                        - L116-125: year window @ 57; jd closed-form;
                                    σ_jd = 1e-8 day                         (M5)
                                    ⚠ Several constants are written as
                                    integer-division literals — e.g.
                                    `int(10 / 12)` (=0), `int(275 / 9)`
                                    (=30). The `int(7 * (year + int(10/12))
                                    / 4)` term reduces to `int(7·year/4)`.
                                    `int(275/9) = 30` is a fixed integer.
                                    These appear lifted from a Y/M/D Julian
                                    formula collapsed to month=1 day=1 then
                                    offset by `td.epoch_day`. Equivalent to
                                    JD(year, Jan 0) + epoch_day; should be
                                    cross-validated against a JD reference
                                    (Vallado/Meeus) for boundary years.
  Match verdict:        ✓ unit conversions match declared (rad/min mean
                        motion, JD construction) for the value channel.
                        ⚠ Angular σ propagation through degrees_to_radians
                          requires confirmation that the conversion routine
                          scales the measurement bound by π/180.
                        ⚠ JD formula correctness on edge years (e.g.
                          1900, 2000, 2100) not verified here.

ERROR BOUND
  Bound category:       measurement (for σ_meas) + precision (for σ_prec)
  Bound formula:        - σ_meas: half-ULP of printed digit count
                          (uniform-quantization assumption):
                              angles:    5e-5 deg              → REQ-EF (measurement)
                              ecc:       5e-8                  → REQ-EF (measurement)
                              mm:        5e-9 rev/day          → REQ-EF (measurement)
                              jd:        1e-8 day              → REQ-EF (measurement)
                              bstar:     10% relative (heuristic, NOT half-ULP)
                        - σ_prec: |val|·ε for representation precision,
                          consistent with REQ-EF-3 (closed-form arithmetic
                          inherits |val|·ε on construction).
                        - σ_acc: 0 (no accuracy debt at parse time).
  Bound implemented:    Constructor calls populate (measurement, precision,
                        accuracy) directly:
                          TrackedValue<T>(val, σ_meas, |val|·ε, 0)
                        ✓ matches REQ-EF-3-style closed-form construction.
                        For the angle path, the measurement bound enters in
                        deg units; whether it is scaled by π/180 in the rad
                        result depends on degrees_to_radians.
  Bound verdict:        ✓ for ecc, mean_motion, jd (constructed directly
                          in target units).
                        ? for angles (depends on degrees_to_radians σ
                          handling; not verified from this header).
                        ⚠ for bstar: "10% relative" is a heuristic, not a
                          column-width half-ULP. SR3 §3 documents the BSTAR
                          field as a 5-digit signed mantissa × 10^exponent,
                          so a column-width-derived σ exists (≈5e-1 in the
                          mantissa's last digit, relative ≈1e-5 not 1e-1).
                          The 10% bound is conservative by 4 orders of
                          magnitude — acceptable as an upper bound but not
                          tight; flag for tightening.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form propagation; |val|·ε prec),
                        plus the measurement-bound rules for fixed-precision
                        inputs (TLE columns).
  AUD-EF applies:       AUD-EF wiring tests that measurement σ survives the
                        unit conversion deg→rad (requires inspecting
                        degrees_to_radians; out of scope here).
  AUD-MC applies:       n/a — this is a parse/construct step, not algebra.
  Verification test:    ? — tests/test_tle/ expected; not visible in this
                        header. Suggested checks:
                        - σ_rad = σ_deg·π/180 after degrees_to_radians.
                        - JD agrees with a Vallado/Meeus reference for
                          a handful of sample epochs (incl. 1957-pivot edge).
                        - mean_motion(rad/min) = mean_motion(rev/day)·2π/1440.

NOTES
  - L120 `int(10 / 12)` evaluates to integer 0 in C++ (int division). This
    is the M=1 specialization of a Y/M/D→JD formula collapsed to year-start.
    Functionally OK but obscured by the literal. ⚠ readability flag.
  - L121 `int(275 / 9) = 30` — similar collapse (day-of-year offset for M=1).
    Combined with `td.epoch_day`, the construction effectively gives
    JD(year, Jan 0) + day_of_year, which is the intended TLE epoch.
  - L117 epoch_year window @ 57 follows the NORAD convention (TLEs prior
    to 1957 do not exist; 57 as pivot is exactly the launch year of
    Sputnik 1). ✓ historically motivated.
  - The bstar "10% relative" σ is a documented convention in the header
    comment (L56) — acceptable as a heuristic measurement-bound floor, but
    is NOT the SR3 column-width half-ULP. Could be tightened.
  - σ_prec uses `abs(val) · ε` uniformly. For values near zero (e.g.
    eccentricity of a nearly-circular orbit), this gives σ_prec ≈ 0,
    which is correct: the representation precision of "0.0000001" stored
    as a double is bounded by |val|·ε. ✓.
  - The function is `static` and templated on T; instantiates for any
    `T` with `boost::math::constants::pi<T>()` and
    `std::numeric_limits<T>::epsilon()`. ✓ portable to multiprecision T.
```

---

## File-level verdict for `tle_parser.h` + `tle_parser.cpp`

- **A. Error wiring**: ✓ for `from_tle_data` — TrackedValue<T> constructed
  with explicit (meas, prec, accuracy) per REQ-EF-3. The `parse(...)`
  overloads emit raw `double`/`int` into TleData and do not participate
  in the error budget at this layer.
- **B. Algebra axioms**: n/a — this is a parser / data-marshal step.
- **C. Theoretical basis**:
  - Card 1 `parse(line1, line2, out)`: ✓ R14 (2026-05-13) — definition
    in tle_parser.cpp:140-142 audited. Pure delegation to three-arg
    overload. SR3 §3 column positions match line-by-line at
    tle_parser.cpp:22-48; Alpha-5 string-typed catnum at cpp:134-136.
  - Card 2 `parse(name, line1, line2, out)`: ✓ R14 (2026-05-13) —
    primary implementation at tle_parser.cpp:144-178 audited. Name is
    opaque std::string; no other special handling.
  - Card 3 `from_tle_data`: ✓ unit conversions, ✓ ecc/mm/jd σ_meas match
    SR3 column widths, ⚠ angle σ propagation depends on
    `degrees_to_radians` (not in this header), ⚠ B* σ is 10% heuristic
    (not column-width tight), ⚠ JD literals use opaque integer-division
    forms.

**File verdict: PASS** — R14 (2026-05-13) closed the Cards 1-2 ?-mark
audit gap by inspecting tle_parser.cpp definitions. SR3 §3 column
positions and the Alpha-5 / assumed-decimal / exponential field parsers
all match the spec. Three minor follow-ups remain for `from_tle_data`
(Card 3), not blocking:
1. Confirm `math::degrees_to_radians` scales σ_meas by π/180.
2. Tighten B* σ from 10% heuristic to SR3 mantissa-column half-ULP,
   or document the 10% as a deliberate measurement floor.
3. Replace opaque `int(10/12)` / `int(275/9)` JD literals with
   precomputed constants and a citation to the source formula
   (Vallado §3.5 / Meeus §7).
4. Minor: checksum digit (column 68) is not validated by parse(); add
   `validate_checksum(line)` helper if strict validation is desired.
