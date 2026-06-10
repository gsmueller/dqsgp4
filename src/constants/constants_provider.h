#pragma once

/// @file constants_provider.h
/// ConstantsProvider: single bundle of fundamental, geodetic, and physical
/// constants injected into the propagator (REQ-CP-1).
///
/// Per CON-1 / AUD-CC-15, the published numerical parameters of each named
/// convention (WGS 84, WGS 72, GRS 80, ...) live in ONE place — the geodesy
/// EquipotentialEllipsoid factories (the reusable single source, C2/F3-a).
/// These provider factories only wrap that ellipsoid in the injection bundle;
/// algorithmic code reads values through the provider, never as literals.
///
/// Sub-bundles (REQ-CP-4..8). v1 carries only:
///
///   earth — geodesy::EquipotentialEllipsoid<T>
///
/// Future phases add `gravity`, `astronomy`, `time`, `fundamentals` as the
/// force lambdas that consume them land.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-9,
///   AUD-CC-10, AUD-CC-12, AUD-CC-15, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../geodesy/equipotential_ellipsoid.h"
#include "../math/tracked_value.h"

namespace constants {

/// Bundle of constants injected into a propagator instance.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct ConstantsProvider {
    geodesy::EquipotentialEllipsoid<T> earth;  ///< Reference ellipsoid plus
                                                ///< derived parameters.

    // --- Constructors ---

    /// Wrap an already-constructed ellipsoid (prefer the named factories below).
    explicit ConstantsProvider(const geodesy::EquipotentialEllipsoid<T>& e)
        : earth(e) {}

    // --- Factory methods ---
    //
    // The defining numbers of each named convention live in ONE place: the
    // geodesy EquipotentialEllipsoid factories (the reusable single source, C2).
    // These provider factories only wrap that ellipsoid in the injection bundle,
    // so WGS 84 / WGS 72 / GRS 80 are never re-typed here (F3-a — the duplication
    // drift risk between this file and the geodesy layer is removed). The
    // resulting ellipsoid is byte-identical to the former inline construction.

    /// WGS 84 (NGA.STND.0036) — delegates to geodesy::EquipotentialEllipsoid::wgs84.
    static ConstantsProvider wgs84(const T& series_tolerance) {
        return ConstantsProvider(
            geodesy::EquipotentialEllipsoid<T>::wgs84(series_tolerance));
    }

    /// WGS 72 (Spacetrack Report #3) — delegates to geodesy::…::wgs72 (J₂-defining).
    static ConstantsProvider wgs72(const T& series_tolerance) {
        return ConstantsProvider(
            geodesy::EquipotentialEllipsoid<T>::wgs72(series_tolerance));
    }

    /// GRS 80 (Moritz 1980) — delegates to geodesy::EquipotentialEllipsoid::grs80.
    static ConstantsProvider grs80(const T& series_tolerance) {
        return ConstantsProvider(
            geodesy::EquipotentialEllipsoid<T>::grs80(series_tolerance));
    }
};

} // namespace constants
