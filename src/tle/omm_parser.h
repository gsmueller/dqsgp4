#pragma once

/// @file omm_parser.h
/// CCSDS Orbit Mean-elements Message (OMM) parser, Key-Value Notation (KVN) —
/// the modern Celestrak / Space-Track element format. It maps an OMM into the
/// SAME `tle::TleData` the two-line parser produces, so `TleElements::
/// from_tle_data` and both propagators consume it unchanged (an alternate
/// front-end, not a second element pipeline).
///
/// EPOCH (ISO-8601) is converted to the TLE epoch_year (2-digit) + epoch_day
/// (fractional day-of-year, leap-aware). Both calendar `YYYY-MM-DDThh:mm:ss` and
/// day-of-year `YYYY-DDDThh:mm:ss` forms are accepted.
///
/// MEAN_MOTION_DOT / _DDOT are taken in the Celestrak GP convention — already
/// the TLE-style n-dot/2 and n-ddot/6 — so they map straight to
/// mean_motion_dt2 / mean_motion_ddt6.
///
/// Header-only and tolerant: a malformed value or a missing mandatory key makes
/// `parse_omm_kvn` return false rather than throw.

#include "tle_parser.h"

#include <map>
#include <sstream>
#include <string>

namespace tle {

namespace omm_detail {

inline std::string trim(const std::string& s) {
    const char* ws = " \t\r\n";
    std::size_t a = s.find_first_not_of(ws);
    if (a == std::string::npos) return "";
    std::size_t b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

/// Parse an OMM ISO epoch into a 2-digit year + fractional day-of-year.
inline bool parse_epoch(const std::string& s, int& year2, double& day_of_year) {
    std::size_t tpos = s.find('T');
    std::string date = (tpos == std::string::npos) ? s : s.substr(0, tpos);
    std::string time = (tpos == std::string::npos) ? "" : s.substr(tpos + 1);

    std::size_t d1 = date.find('-');
    if (d1 == std::string::npos) return false;
    int year = std::stoi(date.substr(0, d1));
    bool leap = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    std::string rest = date.substr(d1 + 1);

    int doy;
    std::size_t d2 = rest.find('-');
    if (d2 == std::string::npos) {
        doy = std::stoi(rest);                       // YYYY-DDD
    } else {
        int month = std::stoi(rest.substr(0, d2));   // YYYY-MM-DD
        int day = std::stoi(rest.substr(d2 + 1));
        if (month < 1 || month > 12) return false;
        static const int cum[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
        doy = cum[month - 1] + day + ((leap && month > 2) ? 1 : 0);
    }

    double frac = 0.0;
    if (!time.empty()) {
        int hh = 0, mm = 0;
        double ss = 0.0;
        std::size_t c1 = time.find(':');
        if (c1 == std::string::npos) {
            hh = std::stoi(time);
        } else {
            hh = std::stoi(time.substr(0, c1));
            std::string r2 = time.substr(c1 + 1);
            std::size_t c2 = r2.find(':');
            if (c2 == std::string::npos) {
                mm = std::stoi(r2);
            } else {
                mm = std::stoi(r2.substr(0, c2));
                ss = std::stod(r2.substr(c2 + 1));
            }
        }
        frac = (hh * 3600.0 + mm * 60.0 + ss) / 86400.0;
    }

    year2 = year % 100;
    day_of_year = static_cast<double>(doy) + frac;
    return true;
}

/// Map an OMM key→value table (from KVN lines OR XML elements) into TleData.
/// Single source for both OMM front-ends (parse_omm_kvn, parse_omm_xml). Returns
/// false on a missing mandatory element or a malformed numeric value.
inline bool populate_from_kv(const std::map<std::string, std::string>& kv, TleData& out) {
    const char* required[] = {"EPOCH",          "MEAN_MOTION",      "ECCENTRICITY",
                              "INCLINATION",     "RA_OF_ASC_NODE",   "ARG_OF_PERICENTER",
                              "MEAN_ANOMALY"};
    for (const char* r : required) {
        if (kv.find(r) == kv.end()) return false;
    }

    auto get = [&](const char* k) -> std::string {
        auto it = kv.find(k);
        return it == kv.end() ? std::string() : it->second;
    };
    auto has = [&](const char* k) { return kv.find(k) != kv.end(); };

    try {
        out.name = get("OBJECT_NAME");
        out.intl_designator = get("OBJECT_ID");
        out.satellite_id = get("NORAD_CAT_ID");
        out.catalog_number = out.satellite_id.empty() ? 0 : decode_alpha5(out.satellite_id);
        out.classification =
            (has("CLASSIFICATION_TYPE") && !get("CLASSIFICATION_TYPE").empty())
                ? get("CLASSIFICATION_TYPE")[0]
                : 'U';

        if (!parse_epoch(get("EPOCH"), out.epoch_year, out.epoch_day)) {
            return false;
        }

        out.mean_motion_rev_day = std::stod(get("MEAN_MOTION"));
        out.eccentricity = std::stod(get("ECCENTRICITY"));
        out.inclination_deg = std::stod(get("INCLINATION"));
        out.raan_deg = std::stod(get("RA_OF_ASC_NODE"));
        out.arg_perigee_deg = std::stod(get("ARG_OF_PERICENTER"));
        out.mean_anomaly_deg = std::stod(get("MEAN_ANOMALY"));

        out.bstar = has("BSTAR") ? std::stod(get("BSTAR")) : 0.0;
        out.mean_motion_dt2 = has("MEAN_MOTION_DOT") ? std::stod(get("MEAN_MOTION_DOT")) : 0.0;
        out.mean_motion_ddt6 = has("MEAN_MOTION_DDOT") ? std::stod(get("MEAN_MOTION_DDOT")) : 0.0;
        out.ephemeris_type = has("EPHEMERIS_TYPE") ? std::stoi(get("EPHEMERIS_TYPE")) : 0;
        out.element_set_number = has("ELEMENT_SET_NO") ? std::stoi(get("ELEMENT_SET_NO")) : 0;
        out.revolution_number = has("REV_AT_EPOCH") ? std::stoi(get("REV_AT_EPOCH")) : 0;
    } catch (...) {
        return false;  // malformed numeric value
    }

    out.line1 = "";
    out.line2 = "";
    out.line1_checksum_valid = true;  // N/A for an OMM
    out.line2_checksum_valid = true;
    return true;
}

} // namespace omm_detail

/// Parse a CCSDS OMM (KVN text) into TleData. Returns false if a mandatory
/// element is missing or a value is malformed. Raw line1/line2 are left empty
/// (an OMM has no two-line text); the checksum flags are set true (N/A).
inline bool parse_omm_kvn(const std::string& text, TleData& out) {
    std::map<std::string, std::string> kv;
    std::istringstream stream(text);
    std::string line;
    while (std::getline(stream, line)) {
        std::size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = omm_detail::trim(line.substr(0, eq));
        std::string val = omm_detail::trim(line.substr(eq + 1));
        if (key.empty() || key == "COMMENT") continue;
        kv[key] = val;
    }

    return omm_detail::populate_from_kv(kv, out);
}

} // namespace tle
