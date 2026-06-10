/// test_omm_kvn (gate E2) — CCSDS OMM (KVN) ingestion into TleData.
///
///   A. every mean-element field maps from its OMM key;
///   B. the ISO EPOCH converts to epoch_year + fractional day-of-year
///      (leap-aware calendar form, and the day-of-year form);
///   C. the parsed OMM feeds TleElements::from_tle_data unchanged (so it drives
///      the same propagators as a TLE);
///   D. a missing mandatory key is rejected.

#include "tle/omm_parser.h"
#include "tle/tle_parser.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using T = double;

namespace {

int failed = 0;

void check(const char* name, bool ok, const std::string& detail = "") {
    if (ok) {
        std::cout << "  PASS: " << name << (detail.empty() ? "" : "  " + detail) << "\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << (detail.empty() ? "" : "  " + detail) << "\n";
    }
}

void check_close(const char* name, double got, double want, double tol) {
    if (std::abs(got - want) <= tol) {
        std::cout << "  PASS: " << name << " = " << std::setprecision(10) << got << "\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << " = " << std::setprecision(17) << got
                  << " (want " << want << ")\n";
    }
}

} // namespace

int main() {
    std::cout << "test_omm_kvn (E2): CCSDS OMM (KVN) ingestion\n\n";

    const std::string omm =
        "CCSDS_OMM_VERS = 3.0\n"
        "COMMENT GENERATED VIA TEST\n"
        "OBJECT_NAME = ISS (ZARYA)\n"
        "OBJECT_ID = 1998-067A\n"
        "CENTER_NAME = EARTH\n"
        "REF_FRAME = TEME\n"
        "TIME_SYSTEM = UTC\n"
        "MEAN_ELEMENT_THEORY = SGP4\n"
        "EPOCH = 2020-12-26T16:00:00.000000\n"
        "MEAN_MOTION = 15.49180000\n"
        "ECCENTRICITY = 0.00012700\n"
        "INCLINATION = 51.6442\n"
        "RA_OF_ASC_NODE = 297.6939\n"
        "ARG_OF_PERICENTER = 119.3897\n"
        "MEAN_ANOMALY = 78.7468\n"
        "EPHEMERIS_TYPE = 0\n"
        "CLASSIFICATION_TYPE = U\n"
        "NORAD_CAT_ID = 25544\n"
        "ELEMENT_SET_NO = 999\n"
        "REV_AT_EPOCH = 25975\n"
        "BSTAR = 0.00012345\n"
        "MEAN_MOTION_DOT = 0.00001234\n"
        "MEAN_MOTION_DDOT = 0.0\n";

    // A. Field mapping.
    std::cout << "=== A. field mapping ===\n";
    tle::TleData td;
    check("parse_omm_kvn ok", tle::parse_omm_kvn(omm, td));
    check("OBJECT_NAME", td.name == "ISS (ZARYA)", td.name);
    check("OBJECT_ID", td.intl_designator == "1998-067A", td.intl_designator);
    check("NORAD_CAT_ID text", td.satellite_id == "25544", td.satellite_id);
    check("catalog_number", td.catalog_number == 25544);
    check("classification", td.classification == 'U');
    check_close("MEAN_MOTION", td.mean_motion_rev_day, 15.4918, 1e-9);
    check_close("ECCENTRICITY", td.eccentricity, 0.000127, 1e-12);
    check_close("INCLINATION", td.inclination_deg, 51.6442, 1e-9);
    check_close("RA_OF_ASC_NODE", td.raan_deg, 297.6939, 1e-9);
    check_close("ARG_OF_PERICENTER", td.arg_perigee_deg, 119.3897, 1e-9);
    check_close("MEAN_ANOMALY", td.mean_anomaly_deg, 78.7468, 1e-9);
    check_close("BSTAR", td.bstar, 0.00012345, 1e-12);
    check_close("MEAN_MOTION_DOT -> dt2", td.mean_motion_dt2, 0.00001234, 1e-12);
    check("ELEMENT_SET_NO", td.element_set_number == 999);
    check("REV_AT_EPOCH", td.revolution_number == 25975);

    // B. EPOCH conversion (Dec 26 2020 is day 361 in a leap year; +16 h).
    std::cout << "\n=== B. EPOCH conversion ===\n";
    check("epoch_year (20)", td.epoch_year == 20);
    check_close("epoch_day (calendar, leap)", td.epoch_day, 361.0 + 16.0 / 24.0, 1e-4);

    const std::string omm_doy =
        "EPOCH = 2021-001T00:00:00\n"
        "MEAN_MOTION = 15.0\nECCENTRICITY = 0.001\nINCLINATION = 51.0\n"
        "RA_OF_ASC_NODE = 0.0\nARG_OF_PERICENTER = 0.0\nMEAN_ANOMALY = 0.0\n";
    tle::TleData td_doy;
    check("parse day-of-year EPOCH", tle::parse_omm_kvn(omm_doy, td_doy));
    check("DOY epoch_year (21)", td_doy.epoch_year == 21);
    check_close("DOY epoch_day (1.0)", td_doy.epoch_day, 1.0, 1e-9);

    // C. The OMM drives the element pipeline.
    std::cout << "\n=== C. feeds from_tle_data ===\n";
    tle::TleElements<T> e = tle::TleElements<T>::from_tle_data(td);
    check("mean_motion finite", std::isfinite(e.mean_motion.value));
    check("inclination finite", std::isfinite(e.inclination.value));
    check("epoch_jd finite", std::isfinite(e.epoch_jd.value));

    // D. Missing mandatory element -> rejected.
    std::cout << "\n=== D. validation ===\n";
    const std::string bad = "MEAN_MOTION = 15.0\nINCLINATION = 51.0\n";  // no EPOCH, etc.
    tle::TleData td_bad;
    check("missing-mandatory rejected", !tle::parse_omm_kvn(bad, td_bad));

    std::cout << "\n" << (failed == 0 ? "PASS" : "FAIL") << " — " << failed
              << " failure(s)\n";
    return failed == 0 ? 0 : 1;
}
