/// test_omm_xml — E3 (DQSGP4 Completion Roadmap, register E3).
///
/// Verifies the CCSDS OMM XML parser (tle/omm_xml_parser.h):
///   - it extracts the mean elements + TLE parameters from the XML;
///   - it produces TleData byte-identical to the KVN parser for equivalent
///     content (the two front-ends share omm_detail::populate_from_kv);
///   - the element-name boundary check correctly distinguishes MEAN_MOTION,
///     MEAN_MOTION_DOT, and MEAN_MOTION_DDOT (the subtle XML scan hazard);
///   - tag attributes (units="…") are skipped;
///   - a missing mandatory element makes it return false (tolerant, no throw);
///   - the result feeds TleElements::from_tle_data.
///
/// ExeGate E3: nonzero exit code on any failed check.

#include "tle/omm_xml_parser.h"
#include "tle/omm_parser.h"
#include "tle/tle_parser.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name << "\n"; }
}

bool close(double a, double b) { return std::abs(a - b) < 1e-9; }

const char* kXml =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<omm id=\"CCSDS_OMM_VERS\" version=\"2.0\">\n"
    " <body><segment>\n"
    "  <metadata>\n"
    "   <OBJECT_NAME>ISS (ZARYA)</OBJECT_NAME>\n"
    "   <OBJECT_ID>1998-067A</OBJECT_ID>\n"
    "  </metadata>\n"
    "  <data>\n"
    "   <meanElements>\n"
    "    <EPOCH>2020-12-26T16:00:00.000</EPOCH>\n"
    "    <MEAN_MOTION units=\"rev/day\">15.49180383</MEAN_MOTION>\n"
    "    <ECCENTRICITY>0.0001724</ECCENTRICITY>\n"
    "    <INCLINATION units=\"deg\">51.6442</INCLINATION>\n"
    "    <RA_OF_ASC_NODE units=\"deg\">339.1770</RA_OF_ASC_NODE>\n"
    "    <ARG_OF_PERICENTER units=\"deg\">69.5021</ARG_OF_PERICENTER>\n"
    "    <MEAN_ANOMALY units=\"deg\">25.0833</MEAN_ANOMALY>\n"
    "   </meanElements>\n"
    "   <tleParameters>\n"
    "    <NORAD_CAT_ID>25544</NORAD_CAT_ID>\n"
    "    <CLASSIFICATION_TYPE>U</CLASSIFICATION_TYPE>\n"
    "    <BSTAR>0.00012345</BSTAR>\n"
    "    <MEAN_MOTION_DOT>0.00001234</MEAN_MOTION_DOT>\n"
    "    <MEAN_MOTION_DDOT>0.0</MEAN_MOTION_DDOT>\n"
    "    <ELEMENT_SET_NO>999</ELEMENT_SET_NO>\n"
    "    <REV_AT_EPOCH>25709</REV_AT_EPOCH>\n"
    "   </tleParameters>\n"
    "  </data>\n"
    " </segment></body>\n"
    "</omm>\n";

// Equivalent KVN of the same OMM (for the XML==KVN equivalence check).
const char* kKvn =
    "OBJECT_NAME = ISS (ZARYA)\n"
    "OBJECT_ID = 1998-067A\n"
    "EPOCH = 2020-12-26T16:00:00.000\n"
    "MEAN_MOTION = 15.49180383\n"
    "ECCENTRICITY = 0.0001724\n"
    "INCLINATION = 51.6442\n"
    "RA_OF_ASC_NODE = 339.1770\n"
    "ARG_OF_PERICENTER = 69.5021\n"
    "MEAN_ANOMALY = 25.0833\n"
    "NORAD_CAT_ID = 25544\n"
    "CLASSIFICATION_TYPE = U\n"
    "BSTAR = 0.00012345\n"
    "MEAN_MOTION_DOT = 0.00001234\n"
    "MEAN_MOTION_DDOT = 0.0\n"
    "ELEMENT_SET_NO = 999\n"
    "REV_AT_EPOCH = 25709\n";

} // namespace

int main() {
    tle::TleData x;
    bool ok = tle::parse_omm_xml(std::string(kXml), x);
    check("OMM XML parses", ok);

    if (ok) {
        // Field mapping.
        check("OBJECT_NAME", x.name == "ISS (ZARYA)");
        check("NORAD_CAT_ID -> satellite_id", x.satellite_id == "25544");
        check("catalog_number decoded", x.catalog_number == 25544);
        check("classification", x.classification == 'U');
        check("mean_motion (attribute skipped)", close(x.mean_motion_rev_day, 15.49180383));
        check("eccentricity", close(x.eccentricity, 0.0001724));
        check("inclination", close(x.inclination_deg, 51.6442));
        check("raan", close(x.raan_deg, 339.1770));
        check("arg_perigee", close(x.arg_perigee_deg, 69.5021));
        check("mean_anomaly", close(x.mean_anomaly_deg, 25.0833));
        check("bstar", close(x.bstar, 0.00012345));

        // The boundary disambiguation: MEAN_MOTION, MEAN_MOTION_DOT, and
        // MEAN_MOTION_DDOT are three distinct elements and must not be confused.
        check("MEAN_MOTION_DOT distinct from MEAN_MOTION", close(x.mean_motion_dt2, 0.00001234));
        check("MEAN_MOTION_DDOT distinct (= 0)", close(x.mean_motion_ddt6, 0.0));

        // 2020 is a leap year; Dec 26 = day 361, +16h = 361.6667 (matches E2).
        check("epoch_year = 20", x.epoch_year == 20);
        check("epoch leap-aware day-of-year ~ 361.667",
              std::abs(x.epoch_day - 361.6666667) < 1e-4);
    }

    // XML and KVN of equivalent content must produce identical TleData.
    tle::TleData k;
    bool kok = tle::parse_omm_kvn(std::string(kKvn), k);
    check("equivalent KVN parses", kok);
    if (ok && kok) {
        check("XML == KVN: satellite_id", x.satellite_id == k.satellite_id);
        check("XML == KVN: epoch_year", x.epoch_year == k.epoch_year);
        check("XML == KVN: epoch_day", x.epoch_day == k.epoch_day);
        check("XML == KVN: mean_motion", x.mean_motion_rev_day == k.mean_motion_rev_day);
        check("XML == KVN: eccentricity", x.eccentricity == k.eccentricity);
        check("XML == KVN: inclination", x.inclination_deg == k.inclination_deg);
        check("XML == KVN: mean_motion_dt2", x.mean_motion_dt2 == k.mean_motion_dt2);
        check("XML == KVN: bstar", x.bstar == k.bstar);
    }

    // Feeds the propagator front-end.
    if (ok) {
        tle::TleElements<double> e = tle::TleElements<double>::from_tle_data(x);
        check("from_tle_data: mean_motion > 0", e.mean_motion.value > 0.0);
        check("from_tle_data: inclination finite", std::isfinite(e.inclination.value));
    }

    // A missing mandatory element -> false (tolerant).
    std::string missing(kXml);
    std::size_t p = missing.find("<EPOCH>");
    std::size_t q = missing.find("</EPOCH>");
    if (p != std::string::npos && q != std::string::npos) {
        missing.erase(p, (q + 8) - p);  // drop the whole EPOCH element
    }
    tle::TleData bad;
    check("missing mandatory EPOCH -> false", !tle::parse_omm_xml(missing, bad));

    std::cout << "\n  OMM XML: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
