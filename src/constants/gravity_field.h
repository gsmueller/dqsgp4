#pragma once

/// @file gravity_field.h
/// One gravity field — the fusion of the zonal column and the tesseral block into
/// a single Stokes-coefficient accessor C(n,m) / S(n,m), as derived in
/// design/derivations/geopotential.md §1. The geopotential (1) has exactly ONE
/// coefficient table; "zonal harmonics" and "tesseral harmonics" are its m = 0
/// column and m ≥ 1 block:
///
///   C(0,0) = 1                      monopole (the defining normalization, exact),
///   C(n,0) = −Jₙ,  S(n,0) = 0       zonal column (Jₙ from the ellipsoid/EGM series),
///   C(n,m) = C_nm, S(n,m) = S_nm    tesseral block (m ≥ 1, denormalized Stokes).
///
/// This class only PRESENTS the existing constants::ZonalHarmonics (as −Jₙ at
/// m = 0) and constants::TesseralHarmonics (m ≥ 1) behind one accessor: every
/// coefficient string and honesty tag is threaded VERBATIM from those two homes
/// (negation preserves the three-error budget), so no provenance categorization
/// moves and CR1B is unaffected. The unified forces/geopotential.h evaluator
/// consumes C/S; zonal-only consumers can still read Jₙ via J(n).

#include "tesseral_harmonics.h"
#include "zonal_harmonics.h"

#include "../math/tracked_value.h"

#include <utility>

namespace constants {

/// A complete gravity field: zonal column ⊕ tesseral block behind C(n,m)/S(n,m).
///
/// @tparam T  Underlying numeric type.
template<typename T>
class GravityField {
public:
    /// Fuse a zonal coefficient set and a tesseral block into one field
    /// (either may be empty; honesty tags pass through verbatim).
    GravityField(ZonalHarmonics<T> zonal, TesseralHarmonics<T> tesseral)
        : zonal_(std::move(zonal)), tesseral_(std::move(tesseral)) {}

    /// Stokes cosine coefficient C_nm: monopole C(0,0)=1, zonal C(n,0)=−Jₙ,
    /// tesseral C(n,m) for m ≥ 1 (zero outside the stored degree/order).
    math::TrackedValue<T> C(int n, int m) const {
        if (m == 0) {
            if (n == 0) return math::exact<T>(1);  // monopole, C_00 ≡ 1 (exact)
            return -zonal_.Jn(n);                  // zonal column C_n0 = −Jₙ
        }
        math::TrackedValue<T> c, s;
        tesseral_.get(n, m, c, s);
        return c;
    }

    /// Stokes sine coefficient S_nm: zero for the m = 0 column, tesseral S(n,m)
    /// for m ≥ 1.
    math::TrackedValue<T> S(int n, int m) const {
        if (m == 0) return math::exact<T>(0);
        math::TrackedValue<T> c, s;
        tesseral_.get(n, m, c, s);
        return s;
    }

    /// Direct zonal coefficient Jₙ (n ≥ 2) for zonal-only consumers.
    math::TrackedValue<T> J(int n) const { return zonal_.Jn(n); }

    /// Highest degree present in either half.
    int max_degree() const {
        const int z = zonal_.max_degree();
        const int t = tesseral_.max_degree();
        return z > t ? z : t;
    }

    /// The zonal half (Jₙ home).
    const ZonalHarmonics<T>& zonal() const { return zonal_; }
    /// The tesseral half (m ≥ 1 block).
    const TesseralHarmonics<T>& tesseral() const { return tesseral_; }

    // --- Standard factories (delegate to the existing coefficient homes) ------

    /// EGM2008 / IERS 2010 — the modern field (zonals + low-degree tesserals).
    static GravityField egm2008(const T& series_tolerance) {
        return GravityField(ZonalHarmonics<T>::egm2008(series_tolerance),
                            TesseralHarmonics<T>::egm2008());
    }

    /// WGS72 / Spacetrack Report #3 — a zonal-only field (no tesserals defined).
    static GravityField wgs72(const T& series_tolerance) {
        return GravityField(ZonalHarmonics<T>::wgs72(series_tolerance),
                            TesseralHarmonics<T>{});
    }

private:
    ZonalHarmonics<T> zonal_;
    TesseralHarmonics<T> tesseral_;
};

}  // namespace constants
