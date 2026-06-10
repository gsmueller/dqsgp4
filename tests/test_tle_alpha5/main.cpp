/// test_tle_alpha5 (gate E1) — Alpha-5 catalog-number decoding and TLE checksum
/// validation, the "updated TLE format" capabilities.
///
///   A. decode_alpha5: classic numbers pass through; Alpha-5 leaders A–Z (with I
///      and O skipped) map to the high part 10–33;
///   B. the col-68 mod-10 checksum validates a correct line and rejects a
///      corrupted one;
///   C. parse() populates catalog_number and the checksum flags for both a
///      classic and an Alpha-5 line.

#include "tle/tle_parser.h"

#include <iostream>
#include <string>

namespace {

int failed = 0;

void check_eq(const char* name, int got, int want) {
    if (got == want) {
        std::cout << "  PASS: " << name << " = " << got << "\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << " = " << got << " (want " << want << ")\n";
    }
}

void check(const char* name, bool ok) {
    if (ok) {
        std::cout << "  PASS: " << name << "\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << "\n";
    }
}

} // namespace

int main() {
    std::cout << "test_tle_alpha5 (E1): Alpha-5 decode + checksum\n\n";

    // A. Alpha-5 decode.
    std::cout << "=== A. decode_alpha5 ===\n";
    check_eq("classic 25544", tle::decode_alpha5("25544"), 25544);
    check_eq("classic 00005", tle::decode_alpha5("00005"), 5);
    check_eq("E8493 (E=14)", tle::decode_alpha5("E8493"), 148493);
    check_eq("A0001 (A=10)", tle::decode_alpha5("A0001"), 100001);
    check_eq("Z9999 (Z=33)", tle::decode_alpha5("Z9999"), 339999);
    check_eq("H0000 (H=17, last before skipped I)", tle::decode_alpha5("H0000"), 170000);
    check_eq("J0000 (J=18, first after skipped I)", tle::decode_alpha5("J0000"), 180000);
    check_eq("N0000 (N=22, last before skipped O)", tle::decode_alpha5("N0000"), 220000);
    check_eq("P0000 (P=23, first after skipped O)", tle::decode_alpha5("P0000"), 230000);
    check_eq("I-leader rejected", tle::decode_alpha5("I1234"), 0);
    check_eq("O-leader rejected", tle::decode_alpha5("O1234"), 0);

    // B. Checksum: construct-verify (robust) + the canonical SGP4-VER sat 5.
    std::cout << "\n=== B. checksum ===\n";
    std::string l1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    std::string l2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

    int c1 = tle::tle_line_checksum(l1);
    std::string l1_fixed = l1.substr(0, 68) + static_cast<char>('0' + c1);
    check("a correctly-stamped line validates", tle::checksum_valid(l1_fixed));

    std::string l1_bad = l1_fixed;
    l1_bad[20] = (l1_bad[20] == '0') ? '1' : '0';  // perturb a body digit
    check("a corrupted line fails the checksum", !tle::checksum_valid(l1_bad));

    std::cout << "  (sat5 L1 computed=" << c1 << " col68='" << l1[68] << "'; L2 computed="
              << tle::tle_line_checksum(l2) << " col68='" << l2[68] << "')\n";
    check("canonical sat5 L1 checksum valid", tle::checksum_valid(l1));
    check("canonical sat5 L2 checksum valid", tle::checksum_valid(l2));

    // C. parse() integration.
    std::cout << "\n=== C. parse() integration ===\n";
    tle::TleData td;
    check("classic parse ok", tle::parse(l1, l2, td));
    check_eq("classic catalog_number", td.catalog_number, 5);
    check("classic checksum flags set",
          td.line1_checksum_valid && td.line2_checksum_valid);

    std::string a1 = l1;
    a1.replace(2, 5, "E8493");  // Alpha-5 catalog in columns 2–6
    tle::TleData ta;
    check("alpha-5 parse ok", tle::parse(a1, l2, ta));
    check_eq("alpha-5 catalog_number", ta.catalog_number, 148493);

    std::cout << "\n" << (failed == 0 ? "PASS" : "FAIL") << " — " << failed
              << " failure(s)\n";
    return failed == 0 ? 0 : 1;
}
