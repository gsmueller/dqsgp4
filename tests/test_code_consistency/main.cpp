/// test_code_consistency — grep-detectable AUD-CC conformance for the
/// dual-quaternion library (AUD-CC-1, -2, -5, -10, -11, -13, -14, -16).
///
/// Scans the five library directories named in design/audit/code_consistency.md
/// — src/math, src/dynamics, src/constants, src/forces, src/integrators — and
/// asserts the mechanically-checkable style/lexicon contract. The remaining
/// AUD-CC items (-3, -4, -6, -8, -9, -12, -15, -17, -18) are "manual review"
/// per the audit's cross-reference table and are not asserted here.
///
/// Exclusion: src/dynamics/state_from_tle.h, state_conversion.h, and
/// dq_sgp4_propagator.h are documented application BRIDGES (they convert a TLE /
/// a StateVector to and from a DQSGP4 State, or wrap the SGP4 seed in a DQ
/// propagator; see design/propagator_choice.md). They are intentionally outside
/// the reusable-library audit — they name SGP4 by design — so skipped.
///
/// Items asserted:
///   AUD-CC-1  every file declares its directory's namespace; no `using
///             namespace`; no anonymous namespace in a header.
///   AUD-CC-2  every header has `#pragma once` and no `#ifndef`/`#define` guard.
///   AUD-CC-5  no LaTeX `$...$` in source (reserved for design/ markdown).
///   AUD-CC-10 `exact_integer` appears only in src/math/tracked_value.h.
///   AUD-CC-11 `boost::math::constants` appears only in angles.h / tracked_value.h.
///   AUD-CC-13 no `using namespace` at any scope (paired with AUD-CC-1).
///   AUD-CC-14 no application coupling (sgp4/tle/brouwer/hoots) in code —
///             scholarly provenance citations in comments are permitted.
///   AUD-CC-16 no line exceeds 100 characters (UTF-8 code points, not bytes).
///
/// Run with the repository root as CWD (the acceptance harness does `& $exe`
/// from the repo root); the source root is located by searching upward from
/// CWD. Exit 0 iff every assertion holds.

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

int checks_passed = 0;
int checks_failed = 0;

void report(const char* aud, const std::vector<std::string>& violations) {
    if (violations.empty()) {
        ++checks_passed;
        std::cout << "  PASS: " << aud << "\n";
    } else {
        ++checks_failed;
        std::cerr << "  FAIL: " << aud << "  (" << violations.size()
                  << " violation(s))\n";
        for (const auto& v : violations) std::cerr << "          " << v << "\n";
    }
}

/// Count Unicode code points in a UTF-8 line: every byte that is not a
/// continuation byte (0b10xxxxxx) begins a new code point. Greek letters and
/// Unicode math symbols in comments are permitted (AUD-CC-5) and must count as
/// one character each toward the 100-column limit (AUD-CC-16).
std::size_t utf8_length(const std::string& s) {
    std::size_t n = 0;
    for (unsigned char c : s) {
        if ((c & 0xC0) != 0x80) ++n;
    }
    return n;
}

std::string lstrip(const std::string& s) {
    std::size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
    return s.substr(i);
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

bool is_ident_char(char c) {
    unsigned char u = static_cast<unsigned char>(c);
    return std::isalnum(u) || c == '_';
}

/// True if `term` occurs in `hay` delimited by non-identifier characters on
/// both sides (so "tle" does not match "subtle" and "state_from_tle").
bool contains_word(const std::string& hay, const std::string& term) {
    std::size_t pos = 0;
    while ((pos = hay.find(term, pos)) != std::string::npos) {
        bool left_ok  = (pos == 0) || !is_ident_char(hay[pos - 1]);
        std::size_t end = pos + term.size();
        bool right_ok = (end >= hay.size()) || !is_ident_char(hay[end]);
        if (left_ok && right_ok) return true;
        pos = end;
    }
    return false;
}

/// Strip a line down to its code part for AUD-CC-14: drop full-line comments
/// (`///`, `//`, or block-comment continuations `*`, `/*`, `*/`) and any
/// trailing `// ...`. Scholarly citations live in comments and are exempt.
std::string code_part(const std::string& line) {
    std::string t = lstrip(line);
    if (t.rfind("///", 0) == 0 || t.rfind("//", 0) == 0 ||
        t.rfind("*", 0) == 0   || t.rfind("/*", 0) == 0 || t.rfind("*/", 0) == 0) {
        return "";
    }
    std::size_t c = line.find("//");
    return (c == std::string::npos) ? line : line.substr(0, c);
}

struct LibDir { std::string rel; std::string ns; };

struct SourceFile {
    fs::path path;
    std::string rel;            // path relative to repo, forward-slashed
    std::string ns;             // expected namespace
    bool is_header = false;
    std::vector<std::string> lines;
};

fs::path find_src_root() {
    const fs::path bases[] = {".", "..", "../..", "../../.."};
    for (const auto& b : bases) {
        if (fs::exists(b / "src" / "math") && fs::exists(b / "src" / "forces")) {
            return b / "src";
        }
    }
    return {};
}

} // anonymous namespace

int main() {
    std::cout << "test_code_consistency: AUD-CC grep-detectable conformance\n\n";

    fs::path src = find_src_root();
    if (src.empty()) {
        std::cerr << "  FAIL: could not locate src/ from CWD "
                  << fs::current_path().string() << "\n";
        return 1;
    }

    const LibDir dirs[] = {
        {"math", "math"}, {"dynamics", "dynamics"}, {"constants", "constants"},
        {"forces", "forces"}, {"integrators", "integrators"}};

    // --- gather in-scope files (skip the documented SGP4 bridge) ---
    std::vector<SourceFile> files;
    for (const auto& d : dirs) {
        fs::path dir = src / d.rel;
        if (!fs::exists(dir)) continue;
        for (const auto& e : fs::recursive_directory_iterator(dir)) {
            if (!e.is_regular_file()) continue;
            const fs::path& p = e.path();
            std::string ext = p.extension().string();
            if (ext != ".h" && ext != ".cpp") continue;
            if (p.filename() == "state_from_tle.h" ||
                p.filename() == "state_conversion.h" ||
                p.filename() == "dq_sgp4_propagator.h") continue;  // application bridges

            SourceFile sf;
            sf.path = p;
            sf.rel = "src/" + d.rel + "/" + p.filename().string();
            sf.ns = d.ns;
            sf.is_header = (ext == ".h");
            std::ifstream in(p);
            std::string line;
            while (std::getline(in, line)) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                sf.lines.push_back(line);
            }
            files.push_back(std::move(sf));
        }
    }
    std::cout << "Scanning " << files.size()
              << " library files (state_from_tle.h, state_conversion.h, "
                 "dq_sgp4_propagator.h excluded)\n\n";

    std::vector<std::string> v_ns, v_using, v_anon, v_pragma, v_guard, v_dollar,
        v_exact, v_boost, v_couple, v_len;
    const std::string couple_terms[] = {"sgp4", "tle", "brouwer", "hoots"};

    for (const auto& f : files) {
        bool has_ns = false, has_pragma = false;
        for (std::size_t i = 0; i < f.lines.size(); ++i) {
            const std::string& line = f.lines[i];
            const std::string stripped = lstrip(line);

            // AUD-CC-16: line length (code points).
            if (utf8_length(line) > 100) {
                v_len.push_back(f.rel + ":" + std::to_string(i + 1) + " (" +
                                std::to_string(utf8_length(line)) + " chars)");
            }
            // AUD-CC-1: namespace presence / AUD-CC-13: using namespace.
            if (line.find("namespace " + f.ns) != std::string::npos) has_ns = true;
            if (line.find("using namespace") != std::string::npos) {
                v_using.push_back(f.rel + ":" + std::to_string(i + 1));
            }
            // AUD-CC-1: no anonymous namespace in a header.
            if (f.is_header && stripped.rfind("namespace {", 0) == 0) {
                v_anon.push_back(f.rel + ":" + std::to_string(i + 1));
            }
            // AUD-CC-2: header guards instead of #pragma once.
            if (stripped == "#pragma once") has_pragma = true;
            if (stripped.rfind("#ifndef", 0) == 0 || stripped.rfind("#define", 0) == 0) {
                // a guard pair is #ifndef X_H / #define X_H; #define of a value
                // is fine, but the library uses #pragma once exclusively, so any
                // #ifndef header guard is a violation.
                if (stripped.rfind("#ifndef", 0) == 0) {
                    v_guard.push_back(f.rel + ":" + std::to_string(i + 1));
                }
            }
            // AUD-CC-5: no LaTeX $...$ in source.
            if (line.find('$') != std::string::npos) {
                v_dollar.push_back(f.rel + ":" + std::to_string(i + 1));
            }
            // AUD-CC-10: exact_integer only in tracked_value.h.
            if (line.find("exact_integer") != std::string::npos &&
                f.path.filename() != "tracked_value.h") {
                v_exact.push_back(f.rel + ":" + std::to_string(i + 1));
            }
            // AUD-CC-11: boost::math::constants only in angles.h / tracked_value.h.
            if (line.find("boost::math::constants") != std::string::npos &&
                f.path.filename() != "angles.h" &&
                f.path.filename() != "tracked_value.h") {
                v_boost.push_back(f.rel + ":" + std::to_string(i + 1));
            }
            // AUD-CC-14: application coupling in code (comments are exempt).
            std::string code = to_lower(code_part(line));
            for (const auto& term : couple_terms) {
                if (contains_word(code, term)) {
                    v_couple.push_back(f.rel + ":" + std::to_string(i + 1) +
                                       " (" + term + ")");
                }
            }
        }
        if (!has_ns) v_ns.push_back(f.rel + " (no `namespace " + f.ns + "`)");
        if (f.is_header && !has_pragma) v_pragma.push_back(f.rel + " (no `#pragma once`)");
    }

    report("AUD-CC-1  namespace declared per file", v_ns);
    report("AUD-CC-1  no anonymous namespace in headers", v_anon);
    report("AUD-CC-1/13 no `using namespace`", v_using);
    report("AUD-CC-2  `#pragma once` in every header", v_pragma);
    report("AUD-CC-2  no `#ifndef` header guards", v_guard);
    report("AUD-CC-5  no LaTeX `$...$` in source", v_dollar);
    report("AUD-CC-10 `exact_integer` only in tracked_value.h", v_exact);
    report("AUD-CC-11 boost constants only in angles.h/tracked_value.h", v_boost);
    report("AUD-CC-14 no application coupling in code", v_couple);
    report("AUD-CC-16 line length <= 100", v_len);

    std::cout << "\n========================================\n";
    std::cout << "Checks passed: " << checks_passed
              << "  failed: " << checks_failed << "\n";
    std::cout << "========================================\n";
    return checks_failed > 0 ? 1 : 0;
}
