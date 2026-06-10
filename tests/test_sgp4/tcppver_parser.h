#pragma once

/// @file tcppver_parser.h
/// @brief Parse the tcppver.out reference output file for SGP4 validation.

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

struct ReferencePoint {
    double t_min;  ///< Time since epoch [minutes]
    double x, y, z;     ///< Position [km]
    double vx, vy, vz;  ///< Velocity [km/s]
};

struct ReferenceCase {
    std::string satellite_id;
    std::vector<ReferencePoint> points;
    /// True when the reference SGP4 errored at t=0 (returned NaN). In tcppver.out
    /// this manifests as a stale buffer: the t=0 short line is byte-identical to
    /// the PREVIOUS satellite's last line (python-sgp4 yields "(Use previous data
    /// line)" in this case — tests.py:727-729). Sat 33334 is the sole such case
    /// in SGP4-VER (perturbed eccentricity out of range, error code 3).
    bool t0_error = false;
};

/// Parse tcppver.out into a vector of reference cases.
/// Format: satellite headers are "NNNNN xx", data rows are space-separated.
inline std::vector<ReferenceCase> parse_tcppver(const std::string& path) {
    std::vector<ReferenceCase> cases;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open " << path << "\n";
        return cases;
    }

    std::string line;
    ReferenceCase current;
    bool have_case = false;
    ReferencePoint prev_last{};      // last point of the previously-completed case
    bool have_prev_last = false;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        // Check for satellite header: ends with " xx"
        if (line.size() >= 3 && line.substr(line.size() - 2) == "xx") {
            if (have_case && !current.points.empty()) {
                prev_last = current.points.back();
                have_prev_last = true;
                cases.push_back(current);
            }
            current = ReferenceCase();
            // Extract satellite ID (everything before " xx")
            std::string trimmed = line.substr(0, line.size() - 3);
            // Trim leading spaces
            size_t start = trimmed.find_first_not_of(' ');
            if (start != std::string::npos) {
                current.satellite_id = trimmed.substr(start);
            } else {
                current.satellite_id = trimmed;
            }
            have_case = true;
            continue;
        }

        if (!have_case) continue;

        // Parse data row: t x y z vx vy vz [optional extra columns]
        std::istringstream iss(line);
        ReferencePoint pt;
        if (iss >> pt.t_min >> pt.x >> pt.y >> pt.z >> pt.vx >> pt.vy >> pt.vz) {
            // Stale-buffer / error detection: a case's FIRST point whose 6 state
            // fields exactly equal the previous case's last point is the
            // "(Use previous data line)" marker — the reference errored at t=0.
            // (Same decimal text parsed identically, so exact == is correct.)
            if (current.points.empty() && have_prev_last
                && pt.x == prev_last.x && pt.y == prev_last.y && pt.z == prev_last.z
                && pt.vx == prev_last.vx && pt.vy == prev_last.vy && pt.vz == prev_last.vz) {
                current.t0_error = true;
            }
            current.points.push_back(pt);
        }
    }

    // Don't forget the last case
    if (have_case && !current.points.empty()) {
        cases.push_back(current);
    }

    return cases;
}
