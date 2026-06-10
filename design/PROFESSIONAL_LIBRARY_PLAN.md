# Professional Library — Theory-First Re-Architecture (Approach)

**★ Governing approach** for the user directive (2026-06-06): *"rebuild the functions so they are extensible,
generalizable and reusable … a professional library … understand each modification from a theoretical basis —
absolutely no modify-first-understand-later … exhaustive theory and knowledge, and the functions just fall
out."* Binding methodology in [[feedback_theory_first_library]] and [[feedback_no_perceived_fidelity]]. This
doc records HOW we approach it; the per-module **theory notes** (in `design/derivations/`) are the substance.

## 1. The method (non-negotiable order)

For every module — new or rebuilt — the sequence is:

1. **Theory note (exhaustive) — FIRST, before any code.** Definitions; governing equations; the *conventions*
   pinned explicitly (frame, epoch, time-scale, units, angle ranges); the physical origin of every term; the
   error budget (the three channels); born-digital references. A claim with no stated convention is a defect.
2. **Derive the abstraction.** The function signatures and the type structure *fall out* of the theory — the
   abstraction mirrors the math. If a function doesn't fall cleanly out of the note, the note is incomplete.
3. **Implement** — the obvious consequence of (2). No discovery happens here.
4. **Isolate + verify against an INDEPENDENT oracle** (no perceived fidelity): each callable unit tested
   alone, matched to truth (astropy/JPL positions, ERFA/analytic identities for transforms, reference
   code/published tables otherwise) to a stated tolerance. Gate it.
5. **Commit** the theory note + code + gate together.

"Absolutely no modify-first." Editing code to see what happens is banned; the understanding precedes the edit.

## 2. Why now — symptom vs. cause

Much of the codebase is *already* professional and should be preserved as the pattern: generic over the
numeric type `T` (double / `cpp_bin_float_50`); `TrackedValue<T>` three-error provenance; pluggable models and
`ModelSelector` presets; a layered tree (`math/` → `constants/ geodesy/ astronomy/` → `forces/ perturbation/`
→ `sgp4/ orbit/ dynamics/`); gated tests.

What needs the treatment is exactly where **the theory was never written down, so conventions leaked** — the
"modify-first" artifacts. The ephemeris we just examined is the textbook symptom: SR3 1970s mean elements on a
1900-based epoch, *no* frame/time-scale abstraction, unwired, and its asserted ±2.5° accuracy never checked
against truth. Nobody could say which frame/epoch `compute_lunar_position` even returns. That is a *cause*
(missing theory), not a bug to patch. The cure is to write the theory and let a correct structure fall out.

## 3. The architecture the theory reveals (worked on the ephemeris)

Writing the geocentric-ephemeris theory exhaustively forces these abstractions into the open — they *fall out*:

```
  r_body(t, target_frame)
      = Frame(source→target, t) · [ Distance(t) · Direction(λ(t), β(t)) ]
  λ,β = mean_longitude(t) + Σ_k  A_k · trig(  i_k·l + j_k·l' + k_k·D + m_k·F + …  )   (Delaunay arguments)
  mean element(t) = element(epoch) + rate·(t − epoch)        [t, epoch in a NAMED time scale]
  A_k, mean elements = born-digital series (VSOP/ELP/Meeus), tagged with provenance + truncation accuracy
```

The reusable, generalizable units this implies:

- **`TimeScale` / `Epoch`** — TT/TDB/UT1/UTC and JD as first-class types, with conversions. (The root fix for
  the epoch confusion: an epoch is never a bare number.)
- **`Frame` + transforms** — ICRF/J2000 ↔ mean-ecliptic-of-date ↔ equatorial, via precession / nutation /
  obliquity (the obliquity generator already exists, SC1). Reusable everywhere; oracle = ERFA + analytic
  identities (orthonormality, round-trip, known angles).
- **`OrbitalModel<Body>`** = {mean elements + rates @ epoch, in a named frame & scale} + {periodic series}
  → (λ, β, r). **Generic over the body** (Sun/Moon/planet share one structure) and **over the element
  SOURCE** (SR3 vs VSOP/ELP vs Meeus) — pluggable exactly like `ModelSelector`. This is what makes O1
  *principled*: the SR3 set and a modern J2000 set are two instances of one theory, each independently
  verifiable, never silently mixed.
- **`Ephemeris`** = `OrbitalModel` + `Frame` + `TimeScale` → a body's position in any frame at any time —
  verified against astropy.

All `T`-generic and `TrackedValue`-tracked. The same Time/Frame/series machinery then serves SGP4, the DQ
propagator, the constants, and any future consumer — i.e. reusable beyond this project.

## 4. Isolation + oracle map (invariant: independent truth, never self-consistency)

| Unit | Independent oracle |
|---|---|
| Time-scale conversions | known epochs / IERS constants; round-trip identities |
| Frame transforms (precession/nutation/obliquity) | ERFA (IAU SOFA) + analytic identities (orthonormal, inverse, J2000 fixed points) |
| Mean-element propagation | published mean longitudes at known dates |
| Periodic series | term-by-term vs the born-digital source; truncation bound vs the next term |
| Full body position | **astropy** geometric ecliptic λ/β,r (and JPL DE if a file is provided) |
| (unchanged) SGP4 path | OR1 frozen oracle — must stay bit-exact |

## 5. Scope & sequence — WHOLE CODEBASE, phased (user, 2026-06-06)

**Scope decided (D-Lib): the entire library, re-architected theory-first, module by module.** Done one layer
at a time, bottom-up by dependency, never modifying before that layer's theory note exists. The authentic
SGP4 path is **OR1-frozen** — it is theory-DOCUMENTED (the BH61/SR3 derivation textbooks already do this), not
rebuilt; only the DQ/library side is restructured, and every step holds 67/67 + OR1 0 km + 33/33.

- **P0 — Foundations audit (math / constants / geodesy).** Largely already professional (`T`-generic,
  `TrackedValue`, honest constants). Action: write/complete their theory notes, confirm each callable is
  isolatable + oracle-gated; rebuild only where a theory note exposes a gap. Mostly ratification.
- **L1 — Time & epoch.** Theory → `TimeScale`/`Epoch` (TT/TDB/TAI/UT1/UTC, JD/MJD, leap seconds, the SGP4 day
  convention). The root of the ephemeris confusion; foundational to everything dated. Oracle: IERS constants,
  round-trip identities, the in-repo SGP4 epoch handling.
- **L2 — Reference frames.** Theory → `Frame` + precession / nutation / obliquity / GMST transforms. Oracle:
  ERFA (IAU SOFA) + analytic identities.
- **L3 — Generic ephemeris.** Theory → `OrbitalModel<Body>` (body-generic, source-pluggable) + `Ephemeris`.
  Oracle: astropy. **Subsumes O1** (SR3 and modern J2000 element sets become two independently-verified
  instances of one theory).
- **L4 — Forces.** Theory → a unified, theory-derived force family: gravity (zonal/tesseral — already strong),
  drag (**subsumes O3 density**, USSA76/NRLMSISE vs the MSIS oracle), third-body (uses L3), SRP. Each
  oracle-verified; the SGP4-path forces stay frozen.
- **L5 — Perturbation / propagation (DQ side).** Theory → the DQ perturbation + integrator wiring as a clean,
  reusable propagator. SGP4 authentic path untouched (OR1).
- **L6 — Integrators & dynamics.** RKF78 / symplectic / RK4 / attitude — theory-document, confirm order
  gates + isolation; rebuild only on a theory-exposed gap.
- **L7 — I/O & API (TLE/OMM, state conversion, ModelSelector).** theory-document + parity gates.

The optional-work items fold in: **O1 → L3, O3 → L4, O2 (EGM Jₙ) → P0** (a constants-provenance theory note).
Existing gates (67/67, OR1, 33/33) are the standing regression backstop throughout.

## 6. Status & replanned sequence (2026-06-09 fresh-eyes replan; supersedes the kickoff "next step")

**Layer status (all gated; `run_acceptance` 73/73, OR1 0 km, test_sgp4 33/33):**

| Layer | Status |
|---|---|
| P0 foundations | partial-implicit: constants covered by the Constants Initiative (8 phases) + formula-layer Stages 0–2; a math/geodesy ratification pass is outstanding (→ R5) |
| L1 time & epoch | ✅ complete (TIME1; `Epoch`/views, value-preserving migration) |
| L2 frames | ✅ complete (FRAME1; rot primitives + sidereal, identity-gated) |
| L3 ephemeris | ✅ complete (EPH + FRAME2; sun/moon Meeus generative series, ecliptic→GCRS chain, JPL-DE430-gated) |
| L4 forces | ~70%: unified geopotential (GEOPOT) + Cartesian third-body force & wiring (TB1) ✅; **real atmosphere (O3) and SRP missing** (→ R1, R2) |
| L5 DQ propagation | facade exists (`DqSgp4Propagator`, G1/F2) but **force injection is missing** — the header's "drag and third-body … layered on via the explicit constructor" is currently FALSE (no force-list parameter; geopotential-only is hardcoded). Presets + adaptive + adapters remain (→ R3) |
| L6 integrators | ✅ complete (RK1: one `ButcherTableau`/`rk_step` family; leapfrog gated + pose-LTE fixed; attitude dispositioned — `runge_kutta_lie_group.md` §7) |
| L7 I/O & API | parsers + parity gates exist (E1/E2/E3, F1/F2, G1); umbrella header, docs/examples, legacy disposition, Earth-fixed (ITRS) outputs remain (→ R4) |

**Replan rationale (the deltas):**
- **Nutation N(t) demoted** from "next" to the R4 Earth-fixed output feature. Its honest forcing consumer is
  an ITRS/ground-track chain; the Moon third-body frame bound it would tighten (~24″ on a ~1e-7 perturbation
  ≈ 1e-11 of total acceleration) is tracked and negligible for propagation. The consumer drives the primitive.
- **"Default-wire third-body?" dissolves into named presets** (R3): the default model stays unchanged
  (OR1-adjacent conservatism); lunisolar/drag/SRP become opt-in presets reachable through the facade —
  which also fixes the doc-lie above.
- **L4 completion promoted to the front**: the largest remaining *honest-fidelity* gaps (the DQ drag path
  carries a 30–50 % Lane/stub accuracy band; SRP is absent), both self-contained with oracles ready — and the
  atmosphere scope is **already user-decided** (OPTIONAL_WORK_PLAN, 2026-06-06: execute O1→O4, staged
  USSA76 → NRLMSISE-00; O1 is done via `moon_meeus`).

**Status (2026-06-09 end of session): R1 ✅ `390ec14` (ATM1) · R2 ✅ `b538b64` (SRP1) — L4 COMPLETE · R3 ✅
`84a861e`+`9e5deeb` (FM1 + AD1/F2) — L5 COMPLETE · R4b ✅ `706a247` (NUT1, the ITRS chain) · R4a ✅ `9565e1d`
(EX1: umbrella `dqsgp4.h`, gated ground-track example, README refresh, legacy disposition, O4 ratified) — L7
COMPLETE · R5 ✅ `0f19335` (P0 ratification note + the O2 EGM2008 Jₙ provenance fix, closes R04; GAL1
double-entry gated). **THE REPLAN IS COMPLETE THROUGH R5 — every roadmap layer P0–L7 is closed.** Remaining:
R6 only (NRLMSISE-00, the flagged multi-session OPT-IN — land on demand; in-repo MSIS_Vers oracle staged,
caller-supplied space weather per decision D2-a). Suite 78/78, OR1 0 km, test_sgp4 33/33.**

**The sequence (each increment: theory note → independent oracle FIRST → impl → standalone-validate → gate → commit):**

- **R1 — USSA76 atmosphere (O3 first stage).** Born-digital 1976 standard-atmosphere layer table → a
  `DensityModel<T>` for the existing `make_drag` seam. Oracle: the published USSA76 table itself
  (table-exact) + in-repo `MSIS_Vers.cpp` sanity points. Replaces the Lane 30–50 % band with a declared
  ~few-% accuracy. OR1-untouched (SGP4 keeps its own Lane drag).
- **R2 — Solar radiation pressure (completes L4).** a = P☉·C_R·(A/m)·(AU/Δ)²·ŝ with a shadow function
  (cylindrical, then conical); **P☉ GENERATED from L☉/(4πc·AU²)** (IAU nominal L☉, defined c and AU — the
  generative-constants mandate); consumes the Sun position (3rd `body_position_gcrs` consumer). Gate:
  generated P☉ vs the published value; shadow analytic identities; GEO/LEO magnitude sanity.
- **R3 — L5 capstone: the facade completes.** `DqSgp4Propagator` force injection + named model presets
  (authentic | boosted | +lunisolar | +drag | +SRP; default unchanged); `Propagator::propagate_adaptive`
  (the standalone RKF78 adaptive loop through the facade — the survey-P5 leftover); the state adapters
  (the km↔m pack written 3× → one module) + `Propagatable` concept. Gates: G1/F2 extensions + a preset gate.
- **R4 — L7 surface.** Umbrella header; README/API refresh; examples (LEO decay with USSA76, GEO lunisolar);
  **legacy-ephemeris disposition** (`solar/lunar_ephemeris.h`, `celestial_body.h` — superseded by the gated
  Meeus instances; document-as-SR3-historical or deprecate, with their gates updated);
  **Earth-fixed output chain**: nutation IAU 2000B (erfa `nut00b` oracle) + equation of equinoxes + polar
  motion + GCRS→ITRS (erfa `c2t06a`) → ground tracks; retro-tightens the third-body frame bound. O4
  (sidereal ratio) folds in.
- **R5 — P0 ratification + O2.** The math/geodesy theory-note ratification pass; O2 = EGM2008 J̄ₙ (n = 5–9)
  provenance re-encode (closes open issue R04).
- **R6 — NRLMSISE-00 (flagged, opt-in, multi-session).** The staged second half of the O3 decision: full port
  vs the in-repo `MSIS_Vers.cpp` oracle, caller-supplied space-weather inputs (decision D2-a). Land on demand.

## 7. The professional-grade close-out (2026-06-09; R6 ✅ `142d4db`+`71f4ed1` — the replan is FULLY complete)

R6 landed as a **mechanical verbatim port** (`gen_msis_port.py`; no hand-typed model arithmetic; MSIS1 measures
the port at machine epsilon vs the reference oracle), with the double-core and single-evaluation-context design
decisions recorded in `nrlmsise00.md` §3. With it, every roadmap layer AND every optional-work item (O1–O4) is
closed; the suite stands at 79 gates.

**The clean-sheet question (user, 2026-06-09: "decide if we need a clean sheet sln file…") — DECIDED: NO
clean-sheet solution; YES generated documentation.**

- A clean-sheet `.sln`/API rewrite is REJECTED. The solution is folder-organized with ~80 projects, every one a
  registered acceptance gate; `run_acceptance.ps1` depends on its layout; the public API already has the layered
  shape a rewrite would aim for (the `dqsgp4.h` umbrella, theory-noted headers, honesty doc-blocks, the
  `Propagatable` concept and presets). Rewriting it would be configuration-control risk (the no-unrecoverable-
  churn rule) for zero functional gain — "immaculate" is achieved by the standing gates, not by a fresh file.
- The genuinely missing professional artifact was the **documentation set** — now generated:
  `tools/gen_docs.py` → `docs/*.html` (16 pages: the overview with the architecture map + the oracle table +
  references; one page per `src/` module with every header's doc-block extracted verbatim and its gates; the
  theory-note index; the acceptance-gate registry parsed from `run_acceptance.ps1`). Generated FROM the
  repository's own headers/notes/gate registry, so the docs regenerate rather than rot.

## 8. NEXT PHASE — documentation excellence + whole-codebase quality pass (user directive, 2026-06-09 close-out)

The user's close-out directive sets the next phase (the Q-phase), in the user's own terms: *"Ensure that the
project enumerates the API as part of the html documentation, with a well integrated help guide. Ensure that the
autogenerated documentation is complete, usable, and up to date, and has textbook quality. Consider having tools
generate integrated diagrams for concepts that would benefit from a diagram. Consider targeted re-writes of any
of the functions, and touch every file; we need to assure the quality of the text and form of each function."*

Decomposed (each increment gated; the standing mandates and frozen invariants apply unchanged):

- **Q1 — API enumeration + help guide.** Extend `gen_docs.py` to parse and enumerate EVERY public
  class/function signature per header (a real API reference, not just doc-blocks), and add an integrated
  task-oriented help guide (propagate a TLE both ways; add forces/presets; adaptive stepping; ITRS/ground
  tracks; reading the three-error budgets; choosing T; space weather) woven into the docs navigation.
  **✅ DONE (2026-06-10).** gen_docs.py gained a STRICT header scanner (an unclassifiable namespace-scope
  construct fails generation): 284 entities / 685 public member signatures across 85 headers → per-header
  API cards on the module pages + the flat `api.html` index (234 public entities; `detail`/`msis`
  enumerated but collapsed as internal). `guide.html` covers the seven mandated tasks with every snippet
  EXTRACTED VERBATIM from gate-compiled sources and every referenced symbol verified against the parsed
  model. Two new gates: **EX2** (`examples/quickstart.cpp` rebuilt as the gated seven-section tour — also
  cures the index claim that called it gate-protected when it had no gate) and **DOC1**
  (`gen_docs.py --check` regenerate-and-diff freshness, negative-tested). Found + recorded en route:
  the DQ integrators' uniform per-step LTE deposit saturates the ACCURACY channel of pose-extracted
  positions (measured: one 30 s RK4 step → 4.3e13 m via r = 2·dual⊗real*; inf over a 60-min arc) —
  documented as the measured semantics in the guide + EX2-pinned; the per-slot redesign is flagged as a
  follow-on task. gen_msis_port.py encoding pinned (the port header is now valid UTF-8; banner-only
  regeneration, MSIS1 green). Suite 79→81 (full), OR1 0 km, 33/33.
- **Q2 — textbook-quality docs.** Every module page carries theory (linked notes with leads), usage (compiling
  snippets), test cases (its gates with what they assert and the measured grades), and references. Docs must
  state MEASURED grades only (the no-perceived-fidelity analog for prose). Completeness check: no public
  symbol undocumented, no stale claim (regenerate-and-diff as the freshness check).
  **✅ DONE (2026-06-10).** Every module page now carries the four sections, each MECHANICAL: Theory = the
  design/derivations notes the module's headers actually cite (title + lead); Usage = snippets extracted
  VERBATIM from the EX1/EX2-gated sources (quickstart grew sections 8–10: constants+Somigliana, Sun/Moon
  GCRS, Kepler/equation-of-centre — including a live tracked-accuracy-MAJORIZES-measured-residual check);
  Test cases = the gates covering the module (doc-block mentions ∪ direct test includes ∪ the I0 foundation
  sweep mirrored from run_acceptance) with their registered measured-grade notes; References = per-module
  citation lists. Doc-completeness RATCHET: 122 undocumented public symbols measured and baselined in
  gen_docs.py (UNDOC_BASELINE; --check fails on regression); Q4's per-module batches lower it to 0, which
  then pins full coverage. Stale-claim fixes en route: equipotential_ellipsoid.h cited nonexistent
  "derivations/001-006" → now cites ch14_equipotential_ellipsoid.md with the real section map. 71/80 fast
  (9 = verifier skips), DOC1 green, OR1 0 km, 33/33.
- **Q3 — generated diagrams.** `tools/gen_diagrams.py` (or gen_docs extension) emitting embedded SVG where a
  diagram genuinely helps: the layer/dependency map, the frame-chain pipeline (TEME/GCRS/ITRS), force
  composition into the propagator, the three-error flow through an operation, the atmosphere model ladder,
  the SGP4-vs-DQSGP4 split. Generated, versioned, regenerable.
  **✅ DONE (2026-06-10).** `tools/gen_diagrams.py` (pure python, no graphviz) emits all six as inline SVG;
  gen_docs embeds them (index: split + layer map; guide: split + frame chain; astronomy/forces/dynamics/
  math/atmosphere module pages) so DOC1's regenerate-and-diff covers diagram freshness automatically. The
  layer map's EDGES are scanned from the real #include graph at generation time (adjacent-layer arrows
  drawn — gutter-confined by construction; the COMPLETE edge list + the four mutual header pairs tabulated
  beside it, all mechanical). Every numeric annotation is a named gate's measured grade. Layout was
  DOM-audited via a live browser (getBBox over every text/box/line): 0 clipped labels, 0 box spills,
  0 arrows through boxes across all six. 71/80 fast (9 = verifier skips), DOC1 green.
- **Q4 — the whole-codebase function-quality pass.** Touch EVERY file: assure the quality of the text and form
  of each function — doc-block accuracy/completeness, naming, comment honesty (constraints stated, no
  narration), dead text, formatting. Targeted REWRITES where form is substandard — value-preserving unless a
  gate covers the change; the OR1-frozen arithmetic is never altered (comment/text improvements in frozen
  files are allowed but verified by the standing bit-gates). Batch per module: fast validation per batch,
  full sweep at boundaries, every gate green throughout.
  **✅ DONE (2026-06-10, two batches).** Batch 1: doc coverage 122→0 undocumented public symbols
  (fix_doc_comments.py promotes //-runs above bare entities — generated files skipped; tracked_value.h's 37
  derivative-bound derivations became real doc-comments; the rest authored; the table generators emit their
  own docs), UNDOC_BASELINE pinned at 0 in DOC1 forever; stale claims fixed (kepler.h's phantom
  equation_convergence_tricks.md → the audited basis note; drag.h "NOT ported"; CustomBuilder's removed drag
  selector; m≠0 "not yet supported" → by-design + the DQ tesseral home); the no-legacy rule made clean across
  src/ (specific implementations named everywhere; verbatim port exempt); ORPHANED src/test_sgp4_ver.cpp
  deleted (zero build references, superseded by tests/test_sgp4). Batch 2: every remaining file read or
  mechanically screened (TODO/stale/narration/dead-code/encoding/dangling-doc-ref sweeps all clean); frames.h's
  "nutation/polar-motion deferred" paragraph corrected to the landed R4b homes; small_angle_series.h's
  per-function docs updated from the old fixed 1e-4 to the B3.2 T-dependent thresholds; gravity_zonal.h's
  "J₃, J₄ … to follow" → the real content + GEOPOT-oracle role.

**★ THE Q-PHASE IS COMPLETE (2026-06-10). FULL SWEEP 81/81 (79 + EX2 + DOC1), OR1 0 km, test_sgp4 33/33
623/623. The documentation set is generated, complete (0 undocumented, freshness-gated), diagrammed, and
guide-integrated; every src/ file passed the quality pass.**
