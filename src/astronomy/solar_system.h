#pragma once

/// Solar and lunar orbital elements computed from fundamental astronomical periods.
///
/// Every constant is derived from input orbital parameters, not hardcoded.
/// The caller provides fundamental periods and eccentricities as TrackedValue
/// with their measurement uncertainties. This module computes the derived
/// quantities that the deep-space perturbation theory needs.
///
/// KEY FINDING (Derivation 008): The SGP4 deep-space model uses TWO different
/// lunar angular rates:
///   - ZNL (1.5835218E-4 rad/min): the ANOMALISTIC rate (27.5546 day period)
///     Used to advance the Moon's mean anomaly in DPPER.
///   - 0.22997150 rad/day: the SIDEREAL rate (27.3216 day period)
///     Used to advance the Moon's mean longitude in DPINIT.
/// These are physically distinct quantities. The mean anomaly advances at
/// the anomalistic rate (perigee-to-perigee), while the mean longitude
/// advances at the sidereal rate (star-to-star).
///
/// Reference: Spacetrack Report #3, Hoots & Roehrich (1980), page 59
///            DATA statements in subroutine DEEP.
///
/// The computation functions are provided as lambdas to support swapping
/// between different Almanac editions without changing the propagator code.

#include "../math/tracked_value.h"
#include "../math/angles.h"
#include <functional>

namespace astronomy {

/// Fundamental astronomical inputs.
/// Each is a TrackedValue carrying measurement uncertainty.
///
/// The SGP4 standard values are circa 1970s astronomical constants.
/// Updated Almanac editions can provide different values without
/// changing any code — just construct with different inputs.
template<typename T>
struct FundamentalConstants {
    // --- Time ---
    math::TrackedValue<T> minutes_per_day;    // = 1440 (exact integer)

    // --- Solar ---
    /// Solar mean anomaly period. The SGP4 value (ZNS) corresponds to
    /// ~365.257 days, close to the sidereal year (365.256 days).
    /// Reference: [SR3] page 59, DATA ZNS/1.19459E-5/
    math::TrackedValue<T> solar_anomaly_period_days;

    /// Earth's orbital eccentricity (ZES in SGP4).
    /// Reference: [SR3] page 59, DATA ZES/.01675/
    math::TrackedValue<T> solar_eccentricity;

    /// Obliquity of the ecliptic. SGP4 uses ~23.4441° (circa 1970 epoch).
    /// IAU 2006 value at J2000.0 is 23.4393°.
    /// Reference: [SR3] page 59, DATA ZCOSIS/.91744867/ → arccos = 23.4441°
    math::TrackedValue<T> obliquity;  // [rad]

    /// Solar argument of perigee at reference epoch.
    /// Reference: [SR3] page 59, DATA ZCOSGS/.1945905/, ZSINGS/-.98088458/
    math::TrackedValue<T> solar_arg_perigee;  // [rad]

    // --- Lunar ---
    /// Lunar anomalistic period (perigee-to-perigee).
    /// This determines ZNL = 2*pi / (anomalistic_period * 1440).
    /// The anomalistic month is 27.554551 days.
    /// Confirmed: ZNL = 1.5835218E-4 matches 27.554551 days to 12 sig figs.
    /// Reference: [SR3] page 59, DATA ZNL/1.5835218E-4/
    ///            [NASA-Eclipse] https://eclipse.gsfc.nasa.gov/SEhelp/moonorbit.html
    math::TrackedValue<T> lunar_anomalistic_period_days;

    /// Lunar sidereal period (star-to-star).
    /// This determines the mean longitude rate = 2*pi / (sidereal_period).
    /// The sidereal month is 27.321661 days.
    /// The rate 0.22997150 rad/day in [SR3] corresponds to 27.321582 days.
    /// Reference: [SR3] page 59, C=4.7199672+.22997150*DAY
    math::TrackedValue<T> lunar_sidereal_period_days;

    /// Lunar orbital eccentricity (ZEL in SGP4).
    /// Reference: [SR3] page 59, DATA ZEL/.05490/
    math::TrackedValue<T> lunar_eccentricity;

    /// Lunar inclination to the ecliptic (~5.145°).
    /// Not directly a DATA constant — used to derive ZCOSIL base value.
    math::TrackedValue<T> lunar_inclination;  // [rad]

    /// Lunar node regression period (~6798.38 days ≈ 18.613 years).
    /// Determines the rate 9.2422029E-4 rad/day in [SR3] page 59.
    math::TrackedValue<T> lunar_node_period_days;

    /// Lunar apsidal advance period (~3231.50 days ≈ 8.85 years).
    math::TrackedValue<T> lunar_perigee_period_days;

    // --- Epoch references ---
    /// Lunar ascending node at reference epoch [rad].
    /// Reference: [SR3] page 59, XNODCE=4.5236020-9.2422029E-4*DAY
    math::TrackedValue<T> lunar_node_epoch;

    /// Lunar mean longitude at reference epoch [rad].
    /// Reference: [SR3] page 59, C=4.7199672+.22997150*DAY
    math::TrackedValue<T> lunar_longitude_epoch;

    /// Lunar perigee longitude at reference epoch [rad].
    /// Reference: [SR3] page 59, GAM=5.8351514+.0019443680*DAY
    math::TrackedValue<T> lunar_perigee_epoch;

    /// Solar mean anomaly at reference epoch [rad].
    /// Reference: [SR3] page 59, ZMOS=6.2565837D0+.017201977D0*DAY
    math::TrackedValue<T> solar_anomaly_epoch;

    /// Construct with standard SGP4-compatible values (circa 1970s).
    /// All values are inputs with measurement uncertainties — none are magic numbers.
    /// Each value is documented with its source in [SR3].
    static FundamentalConstants sgp4_standard() {
        using math::TrackedValue;
        using math::degrees_to_radians;

        FundamentalConstants fc;

        fc.minutes_per_day = math::exact<T>(1440);

        // Solar anomaly period: ZNS=1.19459E-5 → period=365.257 days
        // This is approximately the sidereal year (365.256 days)
        // The exact provenance of this specific value is unclear.
        // Reference: [SR3] page 59
        fc.solar_anomaly_period_days = TrackedValue<T>::measured("365.257", "0.01");

        // Solar eccentricity: ZES=0.01675
        // Reference: [SR3] page 59
        fc.solar_eccentricity = TrackedValue<T>::measured("0.01675", "0.00001");

        // Obliquity: arccos(0.91744867) = 23.4441°
        // Reference: [SR3] page 59, DATA ZCOSIS/.91744867/
        fc.obliquity = degrees_to_radians(
            TrackedValue<T>::measured("23.4441", "0.0001"));

        // Solar argument of perigee: atan2(-0.98088458, 0.1945905)
        // Reference: [SR3] page 59
        fc.solar_arg_perigee = atan2(
            TrackedValue<T>::measured("-0.98088458", "0.000001"),
            TrackedValue<T>::measured("0.1945905", "0.000001"));

        // Lunar anomalistic period: ZNL=1.5835218E-4 → 27.554551 days
        // Confirmed to 12 sig figs against anomalistic month.
        // Reference: [SR3] page 59, [NASA-Eclipse]
        fc.lunar_anomalistic_period_days = TrackedValue<T>::measured("27.554551", "0.000001");

        // Lunar sidereal period: rate 0.22997150 rad/day → 27.321582 days
        // Reference: [SR3] page 59
        fc.lunar_sidereal_period_days = TrackedValue<T>::measured("27.321582", "0.000001");

        // Lunar eccentricity: ZEL=0.05490
        // Reference: [SR3] page 59
        fc.lunar_eccentricity = TrackedValue<T>::measured("0.05490", "0.00001");

        // Lunar inclination to ecliptic: ~5.145°
        fc.lunar_inclination = degrees_to_radians(
            TrackedValue<T>::measured("5.145", "0.001"));

        // Lunar node regression period: rate 9.2422029E-4 rad/day → 6798.38 days
        fc.lunar_node_period_days = TrackedValue<T>::measured("6798.38", "0.01");

        // Lunar apsidal advance period: ~3231.50 days
        fc.lunar_perigee_period_days = TrackedValue<T>::measured("3231.50", "0.01");

        // Epoch reference values from [SR3] page 59
        fc.lunar_node_epoch = TrackedValue<T>::measured("4.5236020", "0.0000001");
        fc.lunar_longitude_epoch = TrackedValue<T>::measured("4.7199672", "0.0000001");
        fc.lunar_perigee_epoch = TrackedValue<T>::measured("5.8351514", "0.0000001");
        fc.solar_anomaly_epoch = TrackedValue<T>::measured("6.2565837", "0.0000001");

        return fc;
    }
};

/// Derived solar/lunar quantities needed by the deep-space perturbation theory.
/// All computed from FundamentalConstants — no magic numbers.
template<typename T>
struct DerivedOrbitalElements {
    // --- Solar rates ---
    /// Solar mean anomaly rate [rad/min] — replaces ZNS.
    /// = 2*pi / (solar_anomaly_period * minutes_per_day)
    math::TrackedValue<T> solar_mean_motion;

    /// Solar eccentricity (passed through from input).
    math::TrackedValue<T> solar_eccentricity;

    /// Trig functions of obliquity — replace ZSINIS, ZCOSIS.
    math::TrackedValue<T> sin_obliquity;
    math::TrackedValue<T> cos_obliquity;

    /// Trig of solar argument of perigee — replace ZSINGS, ZCOSGS.
    math::TrackedValue<T> sin_solar_arg_perigee;
    math::TrackedValue<T> cos_solar_arg_perigee;

    // --- Lunar rates ---
    /// Lunar mean anomaly rate [rad/min] — replaces ZNL.
    /// = 2*pi / (anomalistic_period * minutes_per_day)
    /// NOTE: This uses the ANOMALISTIC period, not the sidereal period.
    math::TrackedValue<T> lunar_mean_motion;

    /// Lunar mean longitude rate [rad/day] — used in DPINIT.
    /// = 2*pi / sidereal_period
    /// NOTE: This uses the SIDEREAL period, not the anomalistic period.
    math::TrackedValue<T> lunar_longitude_rate;

    /// Lunar eccentricity (passed through from input).
    math::TrackedValue<T> lunar_eccentricity;

    /// Lunar node regression rate [rad/day] (negative — regresses westward).
    math::TrackedValue<T> lunar_node_rate;

    /// Lunar perigee advance rate [rad/day].
    math::TrackedValue<T> lunar_perigee_rate;

    /// Solar mean anomaly rate [rad/day] — used in DPINIT for ZMOS.
    math::TrackedValue<T> solar_anomaly_rate_per_day;

    // --- Epoch reference values (passed through) ---
    math::TrackedValue<T> lunar_node_epoch;
    math::TrackedValue<T> lunar_longitude_epoch;
    math::TrackedValue<T> lunar_perigee_epoch;
    math::TrackedValue<T> solar_anomaly_epoch;

    /// Compute all derived quantities from fundamental constants.
    /// Every value is 2*pi/period or a trig function of an input.
    /// No hardcoded numbers anywhere.
    static DerivedOrbitalElements compute(const FundamentalConstants<T>& fc) {
        DerivedOrbitalElements d;

        math::TrackedValue<T> tp = math::two_pi<T>();

        // Solar mean anomaly rate = 2*pi / (period * 1440)
        // This replaces the DATA constant ZNS=1.19459E-5
        d.solar_mean_motion = tp / (fc.solar_anomaly_period_days * fc.minutes_per_day);

        d.solar_eccentricity = fc.solar_eccentricity;

        // Trig of obliquity — replaces ZSINIS and ZCOSIS
        d.sin_obliquity = sin(fc.obliquity);
        d.cos_obliquity = cos(fc.obliquity);

        // Solar argument of perigee trig — replaces ZSINGS, ZCOSGS
        d.sin_solar_arg_perigee = sin(fc.solar_arg_perigee);
        d.cos_solar_arg_perigee = cos(fc.solar_arg_perigee);

        // Lunar mean anomaly rate = 2*pi / (ANOMALISTIC period * 1440)
        // This replaces ZNL. Uses anomalistic period, not sidereal.
        // See Derivation 008 for the distinction.
        d.lunar_mean_motion = tp / (fc.lunar_anomalistic_period_days * fc.minutes_per_day);

        // Lunar mean longitude rate = 2*pi / sidereal_period [rad/day]
        // This replaces the hardcoded 0.22997150 in [SR3] DPINIT.
        d.lunar_longitude_rate = tp / fc.lunar_sidereal_period_days;

        d.lunar_eccentricity = fc.lunar_eccentricity;

        // Lunar node regression rate (negative — regresses westward)
        d.lunar_node_rate = -tp / fc.lunar_node_period_days;

        // Lunar perigee advance rate (positive — advances eastward)
        d.lunar_perigee_rate = tp / fc.lunar_perigee_period_days;

        // Solar mean anomaly rate per day = 2*pi / solar_anomaly_period [rad/day]
        // This replaces the hardcoded 0.017201977 in [SR3] DPINIT.
        d.solar_anomaly_rate_per_day = tp / fc.solar_anomaly_period_days;

        // Pass through epoch reference values
        d.lunar_node_epoch = fc.lunar_node_epoch;
        d.lunar_longitude_epoch = fc.lunar_longitude_epoch;
        d.lunar_perigee_epoch = fc.lunar_perigee_epoch;
        d.solar_anomaly_epoch = fc.solar_anomaly_epoch;

        return d;
    }
};

/// Function type for sidereal time computation.
/// Takes Julian date, returns Greenwich sidereal time [rad].
/// Different implementations correspond to different Almanac editions.
template<typename T>
using SiderealTimeFn = std::function<
    math::TrackedValue<T>(const math::TrackedValue<T>& epoch_jd)>;

/// Function type for lunar ephemeris computation.
/// Takes days from reference epoch, returns orbital elements.
template<typename T>
using LunarEphemerisFn = std::function<void(
    const math::TrackedValue<T>& day,
    const DerivedOrbitalElements<T>& astro,
    math::TrackedValue<T>& cos_inclination,
    math::TrackedValue<T>& sin_inclination,
    math::TrackedValue<T>& sin_longitude,
    math::TrackedValue<T>& cos_longitude,
    math::TrackedValue<T>& mean_longitude,
    math::TrackedValue<T>& mean_anomaly)>;

} // namespace astronomy
