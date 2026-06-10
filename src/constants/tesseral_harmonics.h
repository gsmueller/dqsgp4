#pragma once

/// @file tesseral_harmonics.h
/// Tesseral / sectoral geopotential coefficients C_nm, S_nm (m ≥ 1) — the
/// longitude-DEPENDENT part of the gravity field, complementing the
/// longitude-independent zonals (constants/zonal_harmonics.h). Issue D2.
///
/// The coefficients are stored UNnormalized (the form the Cunningham /
/// Montenbruck-Gill Cartesian recursion in forces/gravity_tesseral.h consumes),
/// built from the born-digital EGM2008 NORMALIZED coefficients C̄_nm, S̄_nm via
/// the geodesy denormalization
///
///   C_nm = N_nm · C̄_nm,   N_nm = √[ (2n+1)(2−δ_0m)(n−m)! / (n+m)! ],
///
/// with √ TRACKED so the precision scales with T (the same generative idea as
/// the zonal Jₙ = −√(2n+1)·C̄ₙ₀ in GAL1). The factorial ratio is built from
/// exact integers, so only the published C̄/S̄ digits limit accuracy.
///
/// Scope (D2 dial: cap at 4×4, EGM2008, with a clear extension path). The
/// canonical degree-2 SECTORAL term C̄_22/S̄_22 (Earth's equatorial ellipticity
/// — the dominant tesseral) is populated born-digital and verified against its
/// unnormalized value. The table and the recursion support any degree/order;
/// adding a higher coefficient is one `add(n, m, "C̄", "S̄")` line sourced from
/// the EGM2008 file (born-digital), so the extension path is mechanical.

#include "../math/tracked_value.h"

#include <string>
#include <vector>

namespace constants {

/// A sparse set of unnormalized tesseral/sectoral Stokes coefficients
/// (m ≥ 1); the zonal m = 0 column lives in ZonalHarmonics.
template<typename T>
struct TesseralHarmonics {
    /// One (n, m) coefficient pair.
    struct Coeff {
        int n;                    ///< degree
        int m;                    ///< order (≥ 1)
        math::TrackedValue<T> C;  ///< unnormalized C_nm
        math::TrackedValue<T> S;  ///< unnormalized S_nm
    };

    std::vector<Coeff> terms;     ///< the stored coefficients (sparse)
    int max_deg = 0;              ///< highest stored degree
    std::string standard;         ///< provenance label (e.g. "egm2008")

    /// Highest stored degree (0 when empty).
    int max_degree() const { return max_deg; }

    /// Unnormalized (C_nm, S_nm); zero for an absent (n, m) (incl. m = 0, which
    /// is the zonal domain handled by ZonalHarmonics).
    void get(int n, int m, math::TrackedValue<T>& C, math::TrackedValue<T>& S) const {
        for (const Coeff& c : terms) {
            if (c.n == n && c.m == m) { C = c.C; S = c.S; return; }
        }
        C = math::exact<T>(0);
        S = math::exact<T>(0);
    }

    /// EGM2008 low-degree sectoral/tesseral coefficients (born-digital).
    static TesseralHarmonics egm2008() {
        TesseralHarmonics th;
        th.standard = "EGM2008";
        // Degree-2 sectoral (equatorial ellipticity): the dominant tesseral.
        // Normalized C̄_22 = 2.43938e-6, S̄_22 = -1.40027e-6 (EGM2008; consistent
        // with the unnormalized C_22 = 1.5746e-6, S_22 = -0.9039e-6 via N_22 =
        // √(5/12) = 0.645497). Source: EGM2008 / IERS Conventions (2010).
        th.add(2, 2, "2.43938e-6", "-1.40027e-6");
        // Extension path: add(3,1,…), add(3,2,…), … add(4,4,…) from the EGM2008
        // coefficient file to fill the 4×4 block (born-digital, one line each).
        return th;
    }

    /// Geodesy denormalization factor N_nm = √[(2n+1)(2−δ_0m)(n−m)!/(n+m)!],
    /// m ≥ 1 (so 2−δ_0m = 2). The factorial ratio (n−m)!/(n+m)! is the exact
    /// reciprocal product 1/∏_{k=n−m+1}^{n+m} k; √ is the tracked sqrt so the
    /// factor's precision scales with T.
    static math::TrackedValue<T> denorm_factor(int n, int m) {
        math::TrackedValue<T> ratio = math::exact<T>((2 * n + 1) * 2);
        for (int k = n - m + 1; k <= n + m; ++k) {
            ratio = ratio / math::exact<T>(k);
        }
        return sqrt(ratio);
    }

    /// Append a coefficient from its born-digital NORMALIZED C̄/S̄ strings,
    /// denormalizing to the unnormalized C_nm/S_nm the recursion consumes.
    void add(int n, int m, const char* Cbar, const char* Sbar) {
        using TV = math::TrackedValue<T>;
        TV f = denorm_factor(n, m);
        // C̄_nm/S̄_nm are gravity-MODEL coefficients (no published formal σ here):
        // model_coefficient books their finite written digits as model-fidelity
        // accuracy and their binary storage as T-scaling precision (2026-06-05
        // panel ruling — uniform, by-nature categorization).
        terms.push_back(Coeff{n, m,
                              f * TV::model_coefficient(Cbar),
                              f * TV::model_coefficient(Sbar)});
        if (n > max_deg) max_deg = n;
    }
};

} // namespace constants
