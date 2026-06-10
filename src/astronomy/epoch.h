#pragma once

/// @file epoch.h
/// L1 — Time scales & epochs (professional-library re-architecture).
/// THEORY FIRST: design/derivations/time_scales_and_epochs.md. This header is
/// the direct consequence of that note's §7 — nothing here is invented.
///
/// A canonical instant is one point on one time line: a (two-part) Julian Date
/// plus a `TimeScale` tag. The library's three historical epoch conventions —
/// a bare JD, the SR3 J1900 day count (JD − 2415020, deep_space.h), and the
/// J2000 Julian century ((JD − 2451545)/36525, sidereal_time.h) — are NOT three
/// epochs; they are VIEWS of one `Epoch`, provided here as pure functions from a
/// single source of the convention bases.
///
/// Convention bases are oracle-anchored against astropy/ERFA (2026-06-06):
///   JD 2451545.0 TT = 2000-01-01 12:00:00  (J2000.0)
///   JD 2415020.0 TT = 1899-12-31 12:00:00  (the J1900 base; exactly 36525 d
///                                           = 100 Julian years before J2000)
/// They are exact-by-convention integer day counts → `exact<T>` (accuracy 0,
/// precision 0), single-sourced here so no consumer re-types a magic offset.
///
/// Scale conversions (TT = TAI+32.184 s; TDB via its periodic series; UT1 via
/// ΔUT1; UTC leap seconds) and calendar↔JD are specified in the theory note
/// (§1, §7) and are added in later L1 increments — each when a consumer needs
/// it and verified against the astropy oracle. They are intentionally NOT
/// pre-implemented here: no unconsumed scaffolding (issue A2).

#include "../math/tracked_value.h"

namespace astronomy {

/// The physical time scales (theory note §1): TAI; TT = TAI + 32.184 s (exact);
/// TDB ≈ TT + periodic; UT1 = UTC + ΔUT1; UTC = TAI − leap seconds. A bare JD
/// must declare which line it lives on; the SGP4 convention treats a TLE/OMM
/// epoch as UTC ≈ UT1 (ΔUT1 = 0).
enum class TimeScale { TAI, TT, TDB, UT1, UTC };

// --- Convention bases (single source of truth; exact by definition) ----------

/// JD 2451545.0 = J2000.0 (2000-01-01 12:00 TT) — exact by convention.
template<typename T> inline math::TrackedValue<T> jd_epoch_j2000()   { return math::exact<T>(2451545); }
/// JD 2415020.0 = the SR3 J1900 base (1899-12-31 12:00 TT) — exact by convention.
template<typename T> inline math::TrackedValue<T> jd_epoch_j1900()   { return math::exact<T>(2415020); }
/// Days per Julian century (36525) — exact by convention.
template<typename T> inline math::TrackedValue<T> days_per_century() { return math::exact<T>(36525); }

/// A canonical instant: a two-part Julian Date (`jd_day` + `jd_frac`, the
/// Vallado/astropy split that keeps sub-day precision — theory note §3) on a
/// named `TimeScale`.
template<typename T>
struct Epoch {
    math::TrackedValue<T> jd_day;   ///< integer / half-integer part of the JD
    math::TrackedValue<T> jd_frac;  ///< sub-day part, |jd_frac| < 0.5
    TimeScale scale;                ///< which time line this instant lives on

    /// Single-JD wrapper (jd_frac = 0). The value-preserving migration target:
    /// `from_jd(epoch_jd)` reproduces the bare-JD arithmetic bit-for-bit, so the
    /// SGP4 path stays OR1-frozen when it adopts the views below.
    static Epoch from_jd(const math::TrackedValue<T>& jd, TimeScale s = TimeScale::UTC) {
        return Epoch{ jd, math::exact<T>(0), s };
    }
    /// Two-part construction (full sub-day precision for new consumers).
    static Epoch from_jd_two_part(const math::TrackedValue<T>& day,
                                  const math::TrackedValue<T>& frac,
                                  TimeScale s = TimeScale::UTC) {
        return Epoch{ day, frac, s };
    }
    /// The combined Julian Date. (Combine only when the magnitude loss is
    /// acceptable; sub-day work should use jd_frac directly — theory note §3.)
    math::TrackedValue<T> jd() const { return jd_day + jd_frac; }
};

// --- Derived VIEWS (pure functions; theory note §5 reconciliation) ------------
// Each reproduces an existing in-repo computation EXACTLY (.value bit-identical),
// so adopting them is value-preserving (OR1-safe).

// JD-form (the fundamental computation, value-IDENTICAL to the in-repo offsets:
// the single-sourced base replaces the inline magic number, with NO extra
// arithmetic, so adopting it is bit-exact / OR1-safe). The Epoch-form delegates.

/// Days from the J1900 base (JD − 2415020) — the SR3 lunisolar argument
/// (replaces deep_space.h's inline `elements.epoch_jd - exact<T>(2415020)`).
template<typename T>
math::TrackedValue<T> days_since_1900(const math::TrackedValue<T>& jd) {
    return jd - jd_epoch_j1900<T>();
}
/// Julian centuries from J2000 ((JD − 2451545)/36525) — the GMST / IAU-polynomial
/// argument (replaces sidereal_time.h's inline century computation).
template<typename T>
math::TrackedValue<T> centuries_since_j2000(const math::TrackedValue<T>& jd) {
    return (jd - jd_epoch_j2000<T>()) / days_per_century<T>();
}

/// Epoch-form: the same view on a typed instant (delegates to the combined JD).
template<typename T>
math::TrackedValue<T> days_since_1900(const Epoch<T>& e) { return days_since_1900(e.jd()); }
/// Epoch-form: the same view on a typed instant (delegates to the combined JD).
template<typename T>
math::TrackedValue<T> centuries_since_j2000(const Epoch<T>& e) { return centuries_since_j2000(e.jd()); }

} // namespace astronomy
