#pragma once

/// @file spherical_harmonics.h
/// The Cunningham (1970) / Montenbruck & Gill real-valued Cartesian recursion for
/// the harmonic functions
///
///   V_nm = (R/r)^{n+1} P_nm(sinφ) cos(mλ),   W_nm = (R/r)^{n+1} P_nm(sinφ) sin(mλ),
///
/// evaluated directly in Earth-fixed Cartesian coordinates — singularity-free (no
/// division by cosφ at the poles) and extensible to any degree/order. This is the
/// pure recursion (no Earth constants), lifted out of forces/gravity_tesseral.h so
/// the zonal (m = 0 column), tesseral (m ≥ 1), and monopole geopotential terms
/// share ONE ladder. Reference: Montenbruck & Gill, "Satellite Orbits" §3.2.4,
/// eqns (3.29)–(3.31). Every step propagates the TrackedValue three-error budget.

#include "tracked_value.h"

#include <vector>

namespace math {

/// V_nm / W_nm tables for 0 ≤ m ≤ n ≤ top, row-major `[n][m]`. The m = 0 column of
/// V is the ZONAL harmonic (R/r)^{n+1} P_n(sinφ); W[n][0] = 0.
template<typename T>
struct VWTable {
    std::vector<std::vector<TrackedValue<T>>> V;
    std::vector<std::vector<TrackedValue<T>>> W;
};

/// Build the V/W tables up to degree `top` from the Earth-fixed position (x,y,z)
/// and Earth radius Re. Bit-identical to the inline recursion it was lifted from.
template<typename T>
VWTable<T> cunningham_vw(const TrackedValue<T>& x, const TrackedValue<T>& y,
                        const TrackedValue<T>& z, const TrackedValue<T>& Re, int top) {
    using TV = TrackedValue<T>;

    TV r2 = x * x + y * y + z * z;
    TV invr2 = exact<T>(1) / r2;
    TV r_mag = sqrt(r2);
    TV xr = x * Re * invr2;   // x·R / r²
    TV yr = y * Re * invr2;   // y·R / r²
    TV zr = z * Re * invr2;   // z·R / r²
    TV RR = Re * Re * invr2;  // R² / r²

    VWTable<T> t;
    t.V.assign(top + 1, std::vector<TV>(top + 1, exact<T>(0)));
    t.W.assign(top + 1, std::vector<TV>(top + 1, exact<T>(0)));
    std::vector<std::vector<TV>>& V = t.V;
    std::vector<std::vector<TV>>& W = t.W;

    V[0][0] = Re / r_mag;  // (R/r), W[0][0] = 0

    // Sectoral diagonal (Montenbruck-Gill 3.29).
    for (int m = 1; m <= top; ++m) {
        TV c = exact<T>(2 * m - 1);
        V[m][m] = c * (xr * V[m - 1][m - 1] - yr * W[m - 1][m - 1]);
        W[m][m] = c * (xr * W[m - 1][m - 1] + yr * V[m - 1][m - 1]);
    }
    // Vertical recursion up each column m (Montenbruck-Gill 3.30/3.31).
    for (int m = 0; m <= top; ++m) {
        for (int n = m + 1; n <= top; ++n) {
            TV a = exact<T>(2 * n - 1) / exact<T>(n - m);
            V[n][m] = a * zr * V[n - 1][m];
            W[n][m] = a * zr * W[n - 1][m];
            if (n - 2 >= m) {
                TV b = exact<T>(n + m - 1) / exact<T>(n - m);
                V[n][m] = V[n][m] - b * RR * V[n - 2][m];
                W[n][m] = W[n][m] - b * RR * W[n - 2][m];
            }
        }
    }
    return t;
}

} // namespace math
