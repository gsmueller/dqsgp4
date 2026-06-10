/// TLE Parser implementation.
///
/// Parses standard two-line element sets into TleData structs.
/// Supports classic and Alpha-5 satellite numbering (decoded to
/// `catalog_number`); the column-68 line checksums are validated and reported.
///
/// Source: Space-Track TLE format specification
///         Vallado (2006) twoline2rv function
///
/// The parser extracts fixed-column fields from the 69-character lines.
/// No magic numbers: every column position and field width is documented.

#include "tle_parser.h"
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace tle {

namespace {

// --- Column positions (0-based) and field widths ---
// Source: Space-Track TLE format, NORAD/USSPACECOM
// Line 1: 69 characters
//   Col  0:    Line number ('1')
//   Col  2-6:  Satellite catalog number (Alpha-5 or numeric)
//   Col  7:    Classification (U/C/S)
//   Col  9-16: International designator
//   Col 18-19: Epoch year (2-digit)
//   Col 20-31: Epoch day of year (fractional)
//   Col 33-42: First derivative of mean motion / 2
//   Col 44-51: Second derivative of mean motion / 6 (decimal point assumed)
//   Col 53-60: BSTAR drag term (decimal point assumed)
//   Col 62:    Ephemeris type
//   Col 64-67: Element set number
//   Col 68:    Checksum
//
// Line 2: 69 characters
//   Col  0:    Line number ('2')
//   Col  2-6:  Satellite catalog number
//   Col  8-15: Inclination [degrees]
//   Col 17-24: RAAN [degrees]
//   Col 26-32: Eccentricity (decimal point assumed)
//   Col 34-41: Argument of perigee [degrees]
//   Col 43-50: Mean anomaly [degrees]
//   Col 52-62: Mean motion [rev/day]
//   Col 63-67: Revolution number at epoch
//   Col 68:    Checksum

// Extract a substring and trim whitespace
std::string extract_field(const std::string& line, int start, int len) {
    if (start + len > static_cast<int>(line.size())) return "";
    std::string s = line.substr(start, len);
    // Trim leading/trailing spaces
    size_t first = s.find_first_not_of(' ');
    if (first == std::string::npos) return "";
    size_t last = s.find_last_not_of(' ');
    return s.substr(first, last - first + 1);
}

// Parse a double from a fixed-width field
double parse_double(const std::string& line, int start, int len) {
    std::string field = extract_field(line, start, len);
    if (field.empty()) return 0.0;
    return std::stod(field);
}

// Parse an integer from a fixed-width field
int parse_int(const std::string& line, int start, int len) {
    std::string field = extract_field(line, start, len);
    if (field.empty()) return 0;
    return std::stoi(field);
}

// Parse the "assumed decimal point" format used for eccentricity
// The TLE stores eccentricity as ".NNNNNNN" without the leading "0."
// e.g., "0012345" means 0.0012345
double parse_assumed_decimal(const std::string& line, int start, int len) {
    std::string field = extract_field(line, start, len);
    if (field.empty()) return 0.0;
    return std::stod("0." + field);
}

// Parse the exponential format used for BSTAR and mean_motion_ddt6
// Format: " NNNNN+E" or " NNNNN-E" where the decimal point is assumed
// e.g., " 12345-4" means 0.12345e-4
double parse_tle_exponential(const std::string& line, int start, int len) {
    std::string field = extract_field(line, start, len);
    if (field.empty()) return 0.0;

    // Find the sign of the exponent (last non-digit character position)
    // Format: [+-]NNNNN[+-]N
    // The mantissa sign may or may not be present
    // The exponent sign is the last +/- before the final digit(s)

    std::string mantissa_str;
    std::string exponent_str;

    // Find the exponent separator (last +/- that isn't the first character)
    int exp_pos = -1;
    for (int i = static_cast<int>(field.size()) - 1; i > 0; --i) {
        if (field[i] == '+' || field[i] == '-') {
            exp_pos = i;
            break;
        }
    }

    if (exp_pos > 0) {
        mantissa_str = field.substr(0, exp_pos);
        exponent_str = field.substr(exp_pos);
    } else {
        mantissa_str = field;
        exponent_str = "0";
    }

    // Insert assumed decimal point if not present
    if (mantissa_str.find('.') == std::string::npos) {
        // Handle leading sign
        if (!mantissa_str.empty() && (mantissa_str[0] == '+' || mantissa_str[0] == '-')) {
            mantissa_str = mantissa_str.substr(0, 1) + "0." + mantissa_str.substr(1);
        } else {
            mantissa_str = "0." + mantissa_str;
        }
    }

    double mantissa = std::stod(mantissa_str);
    int exponent = std::stoi(exponent_str);
    return mantissa * std::pow(10.0, exponent);
}

// Extract the raw catalog-number text field (decoding happens in decode_alpha5).
std::string parse_catalog_number(const std::string& line, int start, int len) {
    return extract_field(line, start, len);
}

} // anonymous namespace

int decode_alpha5(const std::string& field) {
    if (field.empty()) return 0;
    unsigned char c0 = static_cast<unsigned char>(field[0]);
    if (std::isalpha(c0)) {
        // Alpha-5: leading letter A–Z (skip I, O) encodes the high part 10–33.
        char up = static_cast<char>(std::toupper(c0));
        int high;
        if (up >= 'A' && up <= 'H') high = 10 + (up - 'A');
        else if (up >= 'J' && up <= 'N') high = 18 + (up - 'J');
        else if (up >= 'P' && up <= 'Z') high = 23 + (up - 'P');
        else return 0;  // 'I' and 'O' are not valid Alpha-5 leaders
        std::string rest = field.substr(1);
        for (char c : rest) {
            if (c < '0' || c > '9') return 0;
        }
        return rest.empty() ? high * 10000 : high * 10000 + std::stoi(rest);
    }
    for (char c : field) {
        if (c < '0' || c > '9') return 0;
    }
    return std::stoi(field);
}

int tle_line_checksum(const std::string& line) {
    int sum = 0;
    int n = static_cast<int>(line.size());
    if (n > 68) n = 68;  // columns 0..67; the checksum digit at col 68 is excluded
    for (int i = 0; i < n; ++i) {
        char c = line[i];
        if (c >= '0' && c <= '9') sum += c - '0';
        else if (c == '-') sum += 1;
    }
    return sum % 10;
}

bool checksum_valid(const std::string& line) {
    if (line.size() < 69) return false;
    char cs = line[68];
    if (cs < '0' || cs > '9') return false;
    return tle_line_checksum(line) == cs - '0';
}

bool parse(const std::string& line1, const std::string& line2, TleData& out) {
    return parse("", line1, line2, out);
}

bool parse(const std::string& name, const std::string& line1, const std::string& line2, TleData& out) {
    // Validate minimum line lengths
    if (line1.size() < 69 || line2.size() < 69) return false;
    if (line1[0] != '1' || line2[0] != '2') return false;

    out.name = name;
    out.line1 = line1;
    out.line2 = line2;

    // --- Line 1 ---
    out.satellite_id = parse_catalog_number(line1, 2, 5);
    out.catalog_number = decode_alpha5(out.satellite_id);
    out.classification = line1[7];
    out.intl_designator = extract_field(line1, 9, 8);

    out.epoch_year = parse_int(line1, 18, 2);
    out.epoch_day = parse_double(line1, 20, 12);

    out.mean_motion_dt2 = parse_double(line1, 33, 10);
    out.mean_motion_ddt6 = parse_tle_exponential(line1, 44, 8);
    out.bstar = parse_tle_exponential(line1, 53, 8);

    out.ephemeris_type = parse_int(line1, 62, 1);
    out.element_set_number = parse_int(line1, 64, 4);

    // --- Line 2 ---
    out.inclination_deg = parse_double(line2, 8, 8);
    out.raan_deg = parse_double(line2, 17, 8);
    out.eccentricity = parse_assumed_decimal(line2, 26, 7);
    out.arg_perigee_deg = parse_double(line2, 34, 8);
    out.mean_anomaly_deg = parse_double(line2, 43, 8);
    out.mean_motion_rev_day = parse_double(line2, 52, 11);
    out.revolution_number = parse_int(line2, 63, 5);

    // Checksums — validated and reported, not enforced.
    out.line1_checksum_valid = checksum_valid(line1);
    out.line2_checksum_valid = checksum_valid(line2);

    return true;
}

} // namespace tle
