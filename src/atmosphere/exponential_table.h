#pragma once

/**
 * @file exponential_table.h
 * @brief Piecewise-exponential static-atmosphere density — Vallado Table 8-4
 *        (USSA76 / CIRA-72 derived), the table-grade `DensityModel` (replan R1).
 *
 * The published 27-band exponential atmospheric model: within band i,
 *
 *   ρ(h) = ρ₀ᵢ · exp( −(h − h₀ᵢ) / Hᵢ ),     band = largest h₀ᵢ ≤ h,
 *
 * which IS the model's definition (not an approximation of it) — each band's H
 * is the published effective scale height absorbing the in-band temperature
 * profile (hydrostatic equilibrium; theory note §2). Rows transcribed
 * mechanically from the in-repo born-digital source
 * `sgp4_references/vallado_celestrak/software/misc/pascal/ATMOSEXP.DAT`
 * by `tools/gen_atmosphere_table.py` (never typed from memory); the table
 * chains across band edges to max relative residual 9.6e-5 (measured), the
 * transcription-error tripwire the ATM1 gate re-checks.
 *
 * HONESTY (no perceived fidelity — theory note §1): this model claims fidelity
 * to the NAMED PUBLISHED STATIC MODEL (digit floors + the chaining-majorant
 * default band, gate-verified table-exactly) — NOT to the real thermosphere.
 * A static standard cannot represent solar/geomagnetic variability (real LEO
 * density varies ×2–5 over the solar cycle — ERROR SOURCE A-D2 in
 * sgp4_near_earth_drag_theoretical_basis.md §15). That representativeness gap
 * is documented, exposed as the caller-widenable `band_rel_accuracy`, and is
 * what a higher-fidelity space-weather model closes by plugging the same
 * DensityModel interface.
 *
 * Replaces the single-band Lane exponential's declared 30 % (vs USSA76) band on
 * the DQ drag path. The SGP4 path keeps its own frozen Lane power-law drag
 * (atmosphere/density_model.h) — OR1 untouched.
 *
 * Theory: design/derivations/atmosphere_exponential_table.md. Gate: ATM1
 * (tests/test_atmosphere).
 */

#include "../forces/drag.h"          // forces::DensityModel<T> + the make_drag seam
#include "../math/tracked_value.h"

#include <cmath>

namespace atmosphere {

/// One band of the published piecewise-exponential model: base altitude
/// [m, exact integer], nominal density at the base [kg/m³, model_coefficient
/// string], scale height [km, model_coefficient string].
struct ExpAtmosphereRow {
    int h0_m;
    const char* rho0;
    const char* H_km;
};

/// Vallado Table 8-4, transcribed from ATMOSEXP.DAT by gen_atmosphere_table.py.
inline constexpr ExpAtmosphereRow kVallado84[] = {
    {10, "1.225", "7.249"},
    {25000, "3.899e-2", "6.349"},
    {30000, "1.774e-2", "6.682"},
    {40000, "3.972e-3", "7.554"},
    {50000, "1.057e-3", "8.382"},
    {60000, "3.206e-4", "7.714"},
    {70000, "8.770e-5", "6.549"},
    {80000, "1.905e-5", "5.799"},
    {90000, "3.396e-6", "5.382"},
    {100000, "5.297e-7", "5.877"},
    {110000, "9.661e-8", "7.263"},
    {120000, "2.438e-8", "9.473"},
    {130000, "8.484e-9", "12.636"},
    {140000, "3.845e-9", "16.149"},
    {150000, "2.070e-9", "22.523"},
    {180000, "5.464e-10", "29.740"},
    {200000, "2.789e-10", "37.105"},
    {250000, "7.248e-11", "45.546"},
    {300000, "2.418e-11", "53.628"},
    {350000, "9.518e-12", "53.298"},
    {400000, "3.725e-12", "58.515"},
    {450000, "1.585e-12", "60.828"},
    {500000, "6.967e-13", "63.822"},
    {600000, "1.454e-13", "71.835"},
    {700000, "3.614e-14", "88.667"},
    {800000, "1.170e-14", "124.640"},
    {900000, "5.245e-15", "181.045"},
};

/// Number of bands in the published table.
inline constexpr int kVallado84Rows = 27;

/// Tracked density [kg/m³] at geometric altitude `altitude_m` [m above the
/// equatorial radius — the forces/drag.h convention]. Band selection is a value
/// comparison (the density_model.h perigee-regime precedent); below the first
/// base (10 m) the first band extrapolates down, above 900 km the last band
/// extrapolates up (both documented model behavior). `band_rel_accuracy` is the
/// declared band vs the published model — default 1e-3, a 10× majorant of the
/// measured 9.6e-5 chaining residual (theory note §5); widen it for the static-
/// standard-vs-real-thermosphere representativeness (A-D2) as the use demands.
template<typename T>
math::TrackedValue<T> vallado84_density(const math::TrackedValue<T>& altitude_m,
                                        T band_rel_accuracy = T(1) / T(1000)) {
    using TV = math::TrackedValue<T>;
    using std::abs;

    int band = 0;
    for (int i = 1; i < kVallado84Rows; ++i) {
        if (altitude_m.value >= T(kVallado84[i].h0_m)) band = i;
        else break;
    }

    TV h0 = math::exact<T>(kVallado84[band].h0_m);
    TV rho0 = TV::model_coefficient(kVallado84[band].rho0);
    TV H = TV::model_coefficient(kVallado84[band].H_km) * math::exact<T>(1000);

    TV rho = rho0 * exp(-(altitude_m - h0) / H);
    rho = math::add_bound(rho, band_rel_accuracy * abs(rho.value),
                          math::ErrorChannel::accuracy);
    return rho;
}

/// Factory: the table model as a `forces::DensityModel<T>` for `make_drag` —
/// the drop-in upgrade of the single-band Lane exponential on the DQ drag path.
template<typename T>
forces::DensityModel<T> vallado84_density_model(T band_rel_accuracy = T(1) / T(1000)) {
    return [band_rel_accuracy](const math::TrackedValue<T>& altitude_m)
               -> math::TrackedValue<T> {
        return vallado84_density<T>(altitude_m, band_rel_accuracy);
    };
}

} // namespace atmosphere
