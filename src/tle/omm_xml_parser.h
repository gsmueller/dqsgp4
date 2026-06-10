#pragma once

/// @file omm_xml_parser.h
/// CCSDS Orbit Mean-elements Message (OMM) parser, XML serialization — the
/// second CCSDS OMM form alongside the KVN parser in omm_parser.h. It extracts
/// the same named elements from the XML and feeds them through the SAME
/// omm_detail::populate_from_kv mapping, so KVN and XML produce byte-identical
/// `tle::TleData` for equivalent content (one mapping, two front-ends).
///
/// The OMM XML layout is
///   <omm ...><body><segment>
///     <metadata>… OBJECT_NAME, OBJECT_ID …</metadata>
///     <data>
///       <meanElements>EPOCH, MEAN_MOTION, ECCENTRICITY, INCLINATION,
///                     RA_OF_ASC_NODE, ARG_OF_PERICENTER, MEAN_ANOMALY</meanElements>
///       <tleParameters>EPHEMERIS_TYPE, CLASSIFICATION_TYPE, NORAD_CAT_ID,
///                      ELEMENT_SET_NO, REV_AT_EPOCH, BSTAR, MEAN_MOTION_DOT,
///                      MEAN_MOTION_DDOT</tleParameters>
///     </data>
///   </segment></body></omm>
///
/// Header-only and tolerant: a missing mandatory element or a malformed value
/// makes `parse_omm_xml` return false rather than throw. The extractor is a
/// targeted element scan (not a full XML DOM): sufficient for the flat CCSDS OMM
/// schema and dependency-free.

#include "omm_parser.h"

#include <map>
#include <string>

namespace tle {

namespace omm_detail {

/// Extract the text content of the first `<tag>…</tag>` element. The character
/// after the tag name must be a delimiter (`>`, whitespace, or `/`) so a search
/// for `MEAN_MOTION` does NOT match `<MEAN_MOTION_DOT>` — the key disambiguation
/// for the OMM schema. Attributes (`<MEAN_MOTION units="rev/day">`) are skipped.
/// Self-closing/empty elements yield an empty string and are treated as absent.
inline bool extract_element(const std::string& xml, const std::string& tag, std::string& value) {
    const std::string open = "<" + tag;
    std::size_t p = xml.find(open);
    while (p != std::string::npos) {
        std::size_t after = p + open.size();
        char c = (after < xml.size()) ? xml[after] : '\0';
        bool boundary = (c == '>' || c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '/');
        if (boundary) {
            std::size_t gt = xml.find('>', after);
            if (gt == std::string::npos) return false;
            if (gt > 0 && xml[gt - 1] == '/') {       // self-closing: empty
                p = xml.find(open, gt + 1);
                continue;
            }
            const std::string close = "</" + tag + ">";
            std::size_t cpos = xml.find(close, gt + 1);
            if (cpos == std::string::npos) return false;
            value = trim(xml.substr(gt + 1, cpos - (gt + 1)));
            return !value.empty();
        }
        p = xml.find(open, after);
    }
    return false;
}

} // namespace omm_detail

/// Parse a CCSDS OMM (XML text) into TleData via the shared KVN/XML mapping.
/// Returns false on a missing mandatory element or a malformed value.
inline bool parse_omm_xml(const std::string& xml, TleData& out) {
    static const char* tags[] = {
        "OBJECT_NAME", "OBJECT_ID", "NORAD_CAT_ID", "CLASSIFICATION_TYPE",
        "EPOCH", "MEAN_MOTION", "ECCENTRICITY", "INCLINATION", "RA_OF_ASC_NODE",
        "ARG_OF_PERICENTER", "MEAN_ANOMALY", "BSTAR", "MEAN_MOTION_DOT",
        "MEAN_MOTION_DDOT", "EPHEMERIS_TYPE", "ELEMENT_SET_NO", "REV_AT_EPOCH"};

    std::map<std::string, std::string> kv;
    for (const char* t : tags) {
        std::string v;
        if (omm_detail::extract_element(xml, t, v)) {
            kv[t] = v;
        }
    }
    return omm_detail::populate_from_kv(kv, out);
}

} // namespace tle
