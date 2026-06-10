# Specification: Constants Provider

The `constants::ConstantsProvider<T>` interface bundles fundamental,
physical, geodetic, gravitational, astronomical, and time-system constants
that downstream propagator code reads. Swapping providers swaps every
dependent result by construction — no recompile, no global state, no
hidden default (CON-1, CON-2).

Each requirement has a stable identifier (`REQ-CP-N`). Tests live in
`tests/test_constants_provider.cpp`.


## REQ-CP-1 — Single point of injection

Verifies: OBJ-5, CON-2, CON-3.

All Earth-related, gravitational, astronomical, time-system, and
fundamental physical constants are read by propagator code through a
single `ConstantsProvider<T>` instance passed to the propagator at
construction. There is no other path to any of these values from inside
algorithmic code.


## REQ-CP-2 — Factory methods for named conventions

Verifies: OBJ-5, CON-1.

The provider exposes named factory methods that bind the parameters
specified by widely used conventions:

    ConstantsProvider<T>::wgs84(tolerance)
    ConstantsProvider<T>::egm2008(tolerance)
    ConstantsProvider<T>::grs80(tolerance)
    ConstantsProvider<T>::jgm3(tolerance)
    ConstantsProvider<T>::custom(...)

Each factory returns a fully constructed provider. The convention's
published numeric parameters appear only inside the factory body — the
CON-1 magic-number rule is satisfied at the factory boundary.


## REQ-CP-3 — Three-error tracking on every input

Verifies: OBJ-3, CON-5.

Every value stored in a provider is a `TrackedValue<T>` with all three
error categories populated:

- Defined parameters (e.g., WGS 84 a = 6378137 m) have
  measurement = 0 (the value is definitional) and precision > 0 (binary
  representation introduces ULP error).
- Measured constants (e.g., GM, with CODATA-derived uncertainty) have
  measurement > 0 from the authoritative source.
- Accuracy = 0 at provider construction; downstream model truncation
  (J₂-only vs full harmonics, etc.) adds to accuracy as the force model
  is applied.


## REQ-CP-4 — Earth ellipsoid sub-bundle

Verifies: OBJ-5.

The `earth` field is a `geodesy::EquipotentialEllipsoid<T>` constructed
from (a, 1/f, GM, ω, tolerance). The provider does not re-export the
ellipsoid's derived quantities — downstream code reads them through the
ellipsoid's surface (`e2()`, `b()`, `gamma_e()`, `J2n(n)`, …).


## REQ-CP-5 — Gravity field sub-bundle (forthcoming)

Verifies: OBJ-5.

The `gravity` field carries the harmonics that go beyond the ellipsoid:
zonal `{J_n}` for the simple case, full tesseral `{C_nm, S_nm}` for the
EGM-class providers. Force lambdas in `src/forces/` read from this
surface.


## REQ-CP-6 — Astronomy sub-bundle (forthcoming)

Verifies: OBJ-4, OBJ-5.

The `astronomy` field exposes Sun and Moon position lookups for third-body
force lambdas. Interface is a `std::function<Vector3<T>(TrackedValue<T> t)>`
that returns the body's geocentric position at time t.


## REQ-CP-7 — Time system sub-bundle (forthcoming)

Verifies: OBJ-4, OBJ-5.

The `time` field exposes conversions between common time systems (UTC,
TAI, UT1, TT, GPS). Conversions are functions, not stored offsets, since
UT1 in particular requires Earth-orientation parameters that vary in
time.


## REQ-CP-8 — Fundamental constants sub-bundle (forthcoming)

Verifies: OBJ-5.

The `fundamentals` field carries Newton's constant G, the speed of light
c, and other purely physical constants. These do not vary by convention,
but their measured uncertainties are recorded from CODATA.


## REQ-CP-9 — Construction is once-only

Verifies: CON-2, feedback_compute_once.

Provider construction performs all derivations (EquipotentialEllipsoid
constructor work, etc.) once. Accessors are pure reads of stored values.
Propagator lambdas may capture a reference to the provider without
re-derivation cost per call.


## REQ-CP-10 — No global state, no defaults

Verifies: CON-2.

There is no `ConstantsProvider<T>::default_provider()`, no static
default-initialized member, no environment-variable lookup. A propagator
must be given a provider explicitly. Independent propagators in the same
process operate on independent providers.


## Phasing

This specification covers v1 implementation through REQ-CP-4 (Earth
ellipsoid + factory methods for common conventions). REQ-CP-5..8
(gravity, astronomy, time, fundamentals) are specified now but the typed
sub-bundles land as the corresponding force lambdas need them: gravity
first (with the J₂ force), astronomy next (with the third-body force),
time and fundamentals last.


## Cross-reference

Audited by `design/audit/code_consistency.md` (style) and the planned
`tests/test_constants_provider.cpp` (numerical values match published
convention to within stated uncertainty).
