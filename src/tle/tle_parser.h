#pragma once

/// Two-Line Element set parser.
///
/// Parses TLE text lines into numeric orbital elements.
/// Output is TrackedValue<T> with measurement errors derived from
/// the TLE format precision (fixed column widths, known digit counts).
///
/// Supports classic TLE format and Alpha-5 satellite numbering.

#include "../math/tracked_value.h"
#include "../math/angles.h"
#include <string>
#include <functional>

namespace tle {

/// Raw parsed TLE data. All angular values in degrees as read from the TLE.
/// Conversion to radians and TrackedValue happens in to_tracked().
struct TleData {
    std::string name;           // satellite name (line 0, if present)
    std::string line1;          // raw line 1
    std::string line2;          // raw line 2

    // Line 1 fields
    std::string satellite_id;   // catalog number as text (Alpha-5 or numeric)
    int catalog_number;         // decoded numeric catalog (Alpha-5 aware)
    char classification;        // U, C, or S
    std::string intl_designator;
    int epoch_year;             // 2-digit year
    double epoch_day;           // fractional day of year
    double mean_motion_dt2;     // first derivative of mean motion / 2
    double mean_motion_ddt6;    // second derivative of mean motion / 6
    double bstar;               // BSTAR drag term
    int ephemeris_type;         // 0 for SGP4
    int element_set_number;

    // Line 2 fields
    double inclination_deg;
    double raan_deg;
    double eccentricity;        // decimal point assumed (0.xxxxxxx)
    double arg_perigee_deg;
    double mean_anomaly_deg;
    double mean_motion_rev_day; // revolutions per day
    int revolution_number;

    // Validation — computed during parse, NOT enforced (a caller may reject).
    bool line1_checksum_valid;  // line 1 col-68 mod-10 checksum matches
    bool line2_checksum_valid;  // line 2 col-68 mod-10 checksum matches
};

/// Parse TLE lines into TleData.
/// Returns false if parsing fails (malformed lines). The col-68 checksums are
/// validated and reported in the `*_checksum_valid` flags but not enforced.
bool parse(const std::string& line1, const std::string& line2, TleData& out);
bool parse(const std::string& name, const std::string& line1, const std::string& line2, TleData& out);

/// Decode a catalog-number field (classic numeric or Alpha-5) to its integer.
/// Alpha-5: a leading letter A–Z (skipping I and O) encodes the high part 10–33,
/// so "E8493" → 148493; a purely numeric field is parsed as-is ("25544" → 25544).
/// Returns 0 for an empty or malformed field.
int decode_alpha5(const std::string& field);

/// The mod-10 TLE line checksum over columns 0–67 (each digit summed, '-' = 1,
/// every other character 0).
int tle_line_checksum(const std::string& line);

/// True iff the line is ≥ 69 chars and its column-68 digit equals the computed
/// `tle_line_checksum`.
bool checksum_valid(const std::string& line);

/// Convert TleData to TrackedValue orbital elements for use by the propagator.
/// Measurement errors are derived from the TLE format precision:
///   - Angles (i, Ω, ω, M): 4 decimal places in degrees → σ ≈ 0.00005° ≈ 8.7e-7 rad
///   - Eccentricity: 7 decimal digits → σ ≈ 5e-8
///   - Mean motion: 8 decimal digits → σ ≈ 5e-9 rev/day
///   - B*: ~2 significant figures in mantissa → σ ≈ 10% relative
template<typename T>
struct TleElements {
    math::TrackedValue<T> inclination;      // [rad]
    math::TrackedValue<T> raan;             // [rad]
    math::TrackedValue<T> eccentricity;     // dimensionless
    math::TrackedValue<T> arg_perigee;      // [rad]
    math::TrackedValue<T> mean_anomaly;     // [rad]
    math::TrackedValue<T> mean_motion;      // [rad/min]
    math::TrackedValue<T> bstar;            // [1/Earth radii]
    math::TrackedValue<T> epoch_jd;         // Julian date of epoch

    /// Convert from raw TleData to TrackedValue elements.
    static TleElements from_tle_data(const TleData& td) {
        using math::TrackedValue;
        using math::degrees_to_radians;

        TleElements e;

        // ADL-safe abs: works for both double and multiprecision
        using std::abs;

        // Angular measurements: 4 decimal places in degrees
        // Half-digit uncertainty: 0.00005 degrees
        T angle_sigma_deg = T(5) / T(100000);  // 0.00005
        T pi_val = boost::math::constants::pi<T>();
        T eps = std::numeric_limits<T>::epsilon();

        std::function<TrackedValue<T>(double)> make_angle = [&](double deg_val) -> TrackedValue<T> {
            T val(deg_val);
            T prec = abs(val) * eps;
            TrackedValue<T> deg(val, angle_sigma_deg, prec, T(0));
            return degrees_to_radians(deg);
        };

        e.inclination  = make_angle(td.inclination_deg);
        e.raan         = make_angle(td.raan_deg);
        e.arg_perigee  = make_angle(td.arg_perigee_deg);
        e.mean_anomaly = make_angle(td.mean_anomaly_deg);

        // Eccentricity: 7 digits after assumed decimal point
        T ecc_sigma = T(5) / T(100000000); // 5e-8
        T ecc_val(td.eccentricity);
        e.eccentricity = TrackedValue<T>(ecc_val, ecc_sigma, abs(ecc_val) * eps, T(0));

        // Mean motion: convert from rev/day to rad/min
        // Precision: 8 decimal places in rev/day
        T mm_sigma_rev_day = T(5) / T(1000000000); // 5e-9
        T mm_val = T(td.mean_motion_rev_day) * T(2) * pi_val / T(1440);
        T mm_sigma = mm_sigma_rev_day * T(2) * pi_val / T(1440);
        e.mean_motion = TrackedValue<T>(mm_val, mm_sigma, abs(mm_val) * eps, T(0));

        // B* drag term: typically 2 significant figures
        // Relative uncertainty ~10%
        T bstar_val(td.bstar);
        T bstar_sigma = abs(bstar_val) * T(1) / T(10); // 10% relative
        e.bstar = TrackedValue<T>(bstar_val, bstar_sigma, abs(bstar_val) * eps, T(0));

        // Epoch: convert year + fractional day to Julian date
        // The precision of epoch_day (8 decimal places) gives ~0.001 seconds
        int year = td.epoch_year;
        if (year < 57) year += 2000; else year += 1900;
        // Simplified JD computation
        T jd = T(367 * year)
             - T(int(7 * (year + int(10 / 12)) / 4))
             + T(int(275 / 9))
             + T(td.epoch_day)
             + T(1721013) + T(1) / T(2);
        T jd_sigma = T(1) / T(100000000); // 1e-8 days ≈ 0.001 seconds
        e.epoch_jd = TrackedValue<T>(jd, jd_sigma, abs(jd) * eps, T(0));

        return e;
    }
};

} // namespace tle
