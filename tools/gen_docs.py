#!/usr/bin/env python
"""HTML documentation generator (Q-phase: API enumeration + guide, 2026-06-10).

Generates the self-contained static documentation set under docs/ FROM THE
REPOSITORY'S OWN TRUTH, so the docs cannot rot independently of the code:

  - docs/index.html        — overview, architecture map, the oracle table,
                             quick starts, references
  - docs/guide.html        — the task-oriented help guide; every code snippet
                             is extracted VERBATIM from a gate-compiled source
                             (examples/*.cpp), so the guide cannot drift from
                             code that actually builds and runs
  - docs/api.html          — the API index: every namespace-scope entity in
                             every public header, linked to its full entry
  - docs/module_<dir>.html — one page per src/ module: every header with its
                             doc-block, theory-note links, covering gates, and
                             the ENUMERATED API (every public class/function
                             signature, parsed from the header itself)
  - docs/theory.html       — the design/derivations index (title + lead)
  - docs/tests.html        — the acceptance-gate registry parsed from
                             tools/run_acceptance.ps1

Honesty rules (the no-perceived-fidelity analog for documentation):
  - The API reference is PARSED from the headers by a strict scanner: any
    namespace-scope construct the scanner cannot classify is reported and
    FAILS generation (no silent omission).
  - Guide snippets are extracted from `// [guide:tag] ... // [guide:end]`
    regions of gate-compiled sources; a missing marker fails generation.
  - `--check` regenerates everything in memory and diffs against docs/ on
    disk; any drift fails (the DOC1 freshness gate).

Re-run after structural changes: python tools/gen_docs.py
"""
import html
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_diagrams  # noqa: E402  (tool-local import; emits the inline SVGs)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
DERIV = os.path.join(ROOT, "design", "derivations")
DOCS = os.path.join(ROOT, "docs")
EXAMPLES = os.path.join(ROOT, "examples")

MODULES = [
    ("math", "Numeric substrate: TrackedValue three-error framework, vector/quaternion/"
             "dual-quaternion algebra, tracked polynomials and series, Kepler solvers."),
    ("constants", "Honest physical constants: providers, gravity-field coefficients, "
                  "every value provenance-tagged (defined / measured / model_coefficient / generated)."),
    ("geodesy", "The equipotential (level) ellipsoid: defining parameters, derived zonals, "
                "Somigliana normal gravity."),
    ("astronomy", "Time scales and epochs, reference frames, precession and nutation, "
                  "sidereal angles, the GCRS→ITRS chain."),
    ("ephemeris", "Sun and Moon analytical ephemerides (Meeus generative series, DE430-gated) "
                  "and the GCRS position chain; the SR3-historical instances are retained, dispositioned."),
    ("atmosphere", "Atmospheric density: the SGP4 Lane power-law (frozen), the Vallado 8-4 "
                   "static table; a space-weather model plugs the same DensityModel interface."),
    ("forces", "The DQ-side force family: unified geopotential, drag, third-body, solar "
               "radiation pressure — each oracle-gated, each returning a tracked body wrench."),
    ("perturbation", "The analytical perturbation terms of the SGP4 theory (secular "
                     "rates, resonance, Kaula forms, orbit-averaged third-body rates). "
                     "These serve the analytical propagator and are not called directly "
                     "by the numerical one."),
    ("integrators", "The explicit Runge–Kutta family on SE(3): one ButcherTableau + rk_step "
                    "driver (euler/RK4/RKF7(8)) and the symplectic leapfrog."),
    ("dynamics", "The dual-quaternion propagator: state/pose/twist, the force-list engine, "
                 "adaptive stepping, the TLE-seeded DQSGP4 facade with perturbation presets."),
    ("sgp4", "The analytical SGP4/SDP4 propagator, verified against the published "
             "test set (33 satellites, 623 points). Used as the verification reference "
             "for the numerical propagator and for recovering epoch states from TLEs; "
             "include it explicitly (it is not part of the umbrella header)."),
    ("tle", "TLE / OMM (KVN and XML) parsing with checksum and Alpha-5 support."),
    ("orbit", "Orbital-element series helpers (equation of centre, Kepler series)."),
]

# Namespaces that are implementation detail: enumerated, but grouped under a
# collapsed "internal" section on the module page and kept out of api.html.
INTERNAL_NAMESPACES = {"detail", "msis", "omm_detail"}

# Where links OUT of docs/ (the theory notes in design/derivations/) point.
# None  -> relative "../design/derivations/" (a local checkout, the default).
# A URL -> that base instead (a published-site build sets the repository's
#          blob URL so the notes render online; docs/ itself is self-contained
#          either way). The default MUST stay None here: DOC1 pins the
#          committed docs/ to the relative form.
WEB_BASE = "https://github.com/gsmueller/dqsgp4/blob/main/"
DESIGN_BASE = (WEB_BASE + "design/derivations/") if WEB_BASE else "../design/derivations/"

# Doc-completeness ratchet (Q2→Q4): every public symbol carries a doc-comment.
# Reached 0 in the Q4 quality pass (2026-06-10); any new undocumented public
# symbol now fails generation and the DOC1 gate.
UNDOC_BASELINE = 0

# Per-module usage snippets: extracted VERBATIM from gate-compiled sources
# (the [guide:tag] regions), so module usage cannot drift from working code.
QS = "examples/quickstart.cpp"
GT = "examples/ground_track.cpp"
MODULE_USAGE = {
    "math": [(QS, "errors", "EX2"), (QS, "precision", "EX2")],
    "constants": [(QS, "constants", "EX2")],
    "geodesy": [(QS, "geodesy", "EX2")],
    "astronomy": [(GT, "itrs", "EX1")],
    "ephemeris": [(QS, "ephemeris", "EX2")],
    "atmosphere": [],
    "forces": [(QS, "presets", "EX2")],
    "perturbation": [],
    "integrators": [(QS, "integrators", "EX2")],
    "dynamics": [(QS, "facade", "EX2"), (QS, "adaptive", "EX2")],
    "sgp4": [(QS, "sgp4", "EX2")],
    "tle": [(QS, "sgp4", "EX2")],
    "orbit": [(QS, "orbit", "EX2")],
}
MODULE_USAGE_NOTE = {
    "perturbation": "These are the SGP4-internal analytical pieces (secular rates, "
                    "short-period corrections, resonance stepping) — driven by the sgp4 "
                    "module rather than called directly; see the sgp4 module's usage. The "
                    "orbit-averaged third-body instance is SR3-historical and retained "
                    "unwired.",
    "tle": "Parsing is the first step of every propagation; the snippet below (shared "
           "with the sgp4 module) starts from the two TLE lines.",
}

# Per-module reference lists (provenance of the theory and the coefficients —
# citations, not grades; the measured grades live in the gate rows).
# Module-role notes shown at the top of certain module pages.
MODULE_TIER_NOTE = {
    "sgp4": "This module implements the analytical SGP4/SDP4 theory. Within the "
            "library it plays two roles: it is the reference the numerical propagator "
            "is tested against (the test suite reproduces the published verification "
            "set, 33 satellites and 623 points, and pins the results bit-for-bit), "
            "and it recovers epoch states from two-line element sets, whose elements "
            "are defined in terms of this model. It is templated on the numeric type "
            "like the rest of the library. It is not exported by the dqsgp4.h "
            "umbrella header, and an automated check keeps the numerical propagator's "
            "dependence on it confined to the two TLE adapter files.",
    "perturbation": "These headers implement the analytical perturbation terms of "
                    "the SGP4 theory. They are used by the sgp4 module rather than "
                    "called directly by the numerical propagator.",
    "orbit": "Most of this module (element recovery, the modified Kepler solver, "
             "osculating elements, the secular update, and state construction) "
             "belongs to the analytical SGP4 pipeline. kepler_series is the "
             "exception: it provides the Fourier–Bessel anomaly series used by the "
             "ephemeris module.",
}

# Generated diagrams (gen_diagrams.py) embedded per page. The dependency
# map's edges are scanned from the real include graph at generation time.
MODULE_DIAGRAMS = {
    "astronomy": ["frame_chain"],
    "forces": ["force_composition"],
    "dynamics": ["force_composition"],
    "math": ["three_error_flow"],
    "atmosphere": ["atmosphere_ladder"],
}

MODULE_REFS = {
    "math": ["IEEE 754 — the representation-bound semantics of the precision channel.",
             "Boost.Multiprecision (cpp_bin_float_50) and Boost.Math (cyl_bessel_j for the Fourier–Bessel Kepler series).",
             "design/CONSTANTS_INITIATIVE_PLAN.md — the provenance-tagged constant scheme (defined / measured / model_coefficient / generated)."],
    "constants": ["NGA EGM2008 — the in-repo 15-digit normalized-coefficient file (datalib/EGM-08norm100.txt).",
                  "NIMA TR8350.2 — WGS-84 defining parameters; Hoots & Roehrich (1980) — the WGS-72 SGP4 set.",
                  "IAU/IERS adopted values for GM, AU, and the defining rotation constants."],
    "geodesy": ["Moritz (1980), Geodetic Reference System 1980 — the level-ellipsoid closed theory (q0/q0′ series, Somigliana).",
                "NIMA TR8350.2 — WGS-84."],
    "astronomy": ["IERS Conventions (2010), TN36 — precession-nutation, ERA, polar motion.",
                  "IAU SOFA / ERFA — the conformance oracles (pfw06/pmat06, nut06a, gmst06, gst06a, pom00, c2t06a) and the in-repo SOFA nutation table.",
                  "Aoki et al. (1982) — the GMST polynomial of the frozen SGP4 TEME convention.",
                  "Capitaine, Wallace & Chapront (2003) — IAU 2006 precession."],
    "ephemeris": ["Meeus, Astronomical Algorithms (1998) — solar §25 and lunar §47 series (tables generated from pymeeus by tools/gen_lunar_terms.py).",
                  "JPL DE430 — the independent numerical oracle (in-repo Sun/Moon table sunmooneph_430t12.txt).",
                  "Hoots & Roehrich (1980) — the SR3-historical instances (retained, dispositioned)."],
    "atmosphere": ["Vallado, Fundamentals of Astrodynamics — the 8-4 piecewise-exponential table (in-repo born-digital ATMOSEXP.DAT).",
                   "Lane (1965) / Hoots & Roehrich (1980) — the frozen SGP4 power-law density."],
    "forces": ["Montenbruck & Gill, Satellite Orbits — the geopotential acceleration (3.33) and closed-form J₂ (3.34).",
               "Battin, An Introduction to the Mathematics and Methods of Astrodynamics — the cancellation-free f(q) third-body form.",
               "Cunningham (1970) — the V/W recursion.",
               "IAU 2012 (AU) and IAU 2015 B3 (nominal L☉) — the generated SRP pressure constant."],
    "perturbation": ["Hoots & Roehrich (1980), Spacetrack Report #3 — the SGP4 analytical theory.",
                     "Brouwer (1959) — the secular theory; Kaula — the tesseral expansion forms."],
    "integrators": ["Fehlberg (1968), NASA TR R-287 — the RKF7(8) embedded pair.",
                    "Butcher (2003) — order conditions; design/derivations/runge_kutta_lie_group.md — the SE(3) Munthe-Kaas adaptation and the LTE-proxy semantics."],
    "dynamics": ["design/derivations/dq_propagator_facade.md — the facade, presets, and force-injection design.",
                 "Hoots & Roehrich (1980) — the authentic SGP4 epoch seeding (a TLE is WGS-72 by definition).",
                 "The dual-quaternion state design notes (design/derivations/, dynamics section)."],
    "sgp4": ["Hoots & Roehrich (1980), Spacetrack Report #3.",
             "Vallado, Crawford, Hujsak & Kelso (2006), AIAA 2006-6753 — Revisiting Spacetrack Report #3 and the SGP4-VER verification suite (33 satellites, 623 points)."],
    "tle": ["The Space-Track TLE format specification (incl. Alpha-5).",
            "CCSDS 502.0-B Orbit Data Messages — the OMM KVN/XML forms."],
    "orbit": ["design/derivations/ephemeris_series.md — the Fourier–Bessel solution of Kepler's equation (converges for all e < 1).",
              "Hoots & Roehrich (1980) — the modified-Kepler iteration (U = E + ω) of the SGP4 path."],
}

CSS = """
body { font-family: Segoe UI, Helvetica, Arial, sans-serif; margin: 0; color: #1a232e;
       background: #f6f8fa; line-height: 1.55; }
header { background: #102a43; color: #f0f4f8; padding: 18px 36px; }
header h1 { margin: 0 0 4px 0; font-size: 1.45em; }
header p { margin: 0; color: #bcccdc; }
nav { background: #243b53; padding: 8px 36px; }
nav a { color: #d9e2ec; margin-right: 18px; text-decoration: none; font-size: 0.95em; }
nav a:hover { color: #fff; text-decoration: underline; }
main { max-width: 1080px; margin: 24px auto; padding: 0 24px; }
h2 { color: #102a43; border-bottom: 2px solid #d9e2ec; padding-bottom: 4px; }
h3 { color: #243b53; margin-bottom: 4px; }
table { border-collapse: collapse; width: 100%; margin: 12px 0; background: #fff; }
th, td { border: 1px solid #d9e2ec; padding: 6px 10px; text-align: left;
         vertical-align: top; font-size: 0.93em; }
th { background: #e3ecf3; }
code, pre { font-family: Consolas, monospace; background: #eef2f6; border-radius: 3px; }
code { padding: 1px 4px; }
pre { padding: 10px 14px; overflow-x: auto; border: 1px solid #d9e2ec; }
.card { background: #fff; border: 1px solid #d9e2ec; border-radius: 6px;
        padding: 12px 18px; margin: 10px 0; }
.muted { color: #627d98; font-size: 0.9em; }
.ent { background: #fff; border: 1px solid #d9e2ec; border-left: 4px solid #486581;
       border-radius: 4px; padding: 8px 14px; margin: 10px 0; }
.ent pre.sig { background: #f0f4f8; border: none; margin: 0 0 6px 0; padding: 8px 10px;
               white-space: pre-wrap; }
.entdoc { margin: 4px 0; }
.kind { display: inline-block; background: #486581; color: #fff; border-radius: 3px;
        font-size: 0.78em; padding: 1px 7px; margin-right: 8px; vertical-align: middle; }
.member-tbl td { font-size: 0.9em; }
details { margin: 10px 0; }
details summary { cursor: pointer; color: #486581; font-weight: 600; }
.provenance { color: #627d98; font-size: 0.85em; margin-top: 4px; }
.docblock { white-space: pre-wrap; border: none; background: none; margin: 0; padding: 0; }
footer { text-align: center; color: #627d98; font-size: 0.85em; margin: 28px 0; }
"""


def page(title, body):
    nav = ('<a href="index.html">Overview</a><a href="guide.html">Guide</a>'
           '<a href="api.html">API index</a><a href="theory.html">Theory notes</a>'
           '<a href="tests.html">Acceptance gates</a>'
           + "".join(f'<a href="module_{m}.html">{m}</a>' for m, _ in MODULES
                     if os.path.isdir(os.path.join(SRC, m))))
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>{html.escape(title)} — DQSGP4</title><style>{CSS}</style></head>
<body><header><h1>{html.escape(title)}</h1>
<p>Arbitrary-precision SGP4 &amp; dual-quaternion orbit propagation — generated documentation</p></header>
<nav>{nav}</nav><main>
{body}
</main><footer>Generated by tools/gen_docs.py from the repository's own headers, theory notes
and gate registry — regenerate after structural changes.</footer></body></html>"""


# ============================ C++ header parsing =============================
#
# A strict line-oriented scanner for the house style: namespace-scope entities
# (struct/class/enum/using/function/variable) with their /// doc-comments, and
# the PUBLIC members of each struct/class. Anything at namespace or class scope
# the scanner cannot classify is recorded as a warning, and warnings fail the
# generation — completeness is enforced, not assumed.

ID_RE = r"[A-Za-z_]\w*"


def _mask_code(line, in_block):
    """Split a raw line into (masked_code, raw_code, doc_comment, in_block).

    masked_code has string/char-literal contents replaced by 'x' (length-
    preserving) and comments removed, so structural scans can count braces and
    parens safely. raw_code is the untouched code portion (same length as
    masked_code). doc_comment is the text of a ``///`` comment if present.
    Block comments are tracked across lines via in_block."""
    masked = []
    raw = []
    doc = None
    i, n = 0, len(line)
    while i < n:
        if in_block:
            j = line.find("*/", i)
            if j < 0:
                return "".join(masked), "".join(raw), doc, True
            i = j + 2
            continue
        c = line[i]
        if c == '"':
            masked.append('"'); raw.append('"'); i += 1
            while i < n and line[i] != '"':
                if line[i] == "\\" and i + 1 < n:
                    masked.append("xx"); raw.append(line[i:i + 2]); i += 2
                else:
                    masked.append("x"); raw.append(line[i]); i += 1
            if i < n:
                masked.append('"'); raw.append('"'); i += 1
        elif c == "'":
            # Char literal (only valid in code); scan for a nearby close.
            j = line.find("'", i + 1)
            if 0 < j <= i + 4:
                masked.append("'" + "x" * (j - i - 1) + "'")
                raw.append(line[i:j + 1])
                i = j + 1
            else:
                masked.append(c); raw.append(c); i += 1
        elif c == "/" and i + 1 < n and line[i + 1] == "/":
            comment = line[i:]
            if comment.startswith("///"):
                doc = comment
            break
        elif c == "/" and i + 1 < n and line[i + 1] == "*":
            j = line.find("*/", i + 2)
            if j < 0:
                return "".join(masked), "".join(raw), doc, True
            i = j + 2
        else:
            masked.append(c); raw.append(c); i += 1
    return "".join(masked), "".join(raw), doc, False


class HeaderParser:
    """Parse one header into a list of entity dicts."""

    def __init__(self, path, rel):
        self.rel = rel
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()
        self.file_doc = self._extract_file_doc(text)
        # Strip /* ... */ blocks up front (newline-preserving) and remember
        # each block's text keyed by its END line, so a block doc-comment
        # above a declaration attaches exactly like a /// run. The first
        # /** block is the @file doc and is not re-attached.
        block_docs = {}
        first_block = True

        def _strip_block(m):
            nonlocal first_block
            body = m.group(0)
            end_line = text.count("\n", 0, m.end())
            if not first_block:
                lines = [re.sub(r"^\s*\*+\s?", "", ln)
                         for ln in body.strip("/*").splitlines()]
                block_docs[end_line] = "\n".join(lines).strip()
            first_block = False
            return re.sub(r"[^\n]", " ", body)

        text = re.sub(r"/\*.*?\*/", _strip_block, text, flags=re.S)
        self.lines = []          # (masked, raw, doc_comment)
        in_block = False
        for idx, ln in enumerate(text.splitlines()):
            masked, raw, doc, in_block = _mask_code(ln, in_block)
            if doc is None and idx in block_docs and block_docs[idx]:
                doc = "///" + block_docs[idx].replace("\n", "\n///")
            self.lines.append((masked, raw, doc))
        self.i = 0
        self.ns = []             # namespace stack
        self.entities = []
        self.warnings = []

    @staticmethod
    def _extract_file_doc(text):
        m = re.search(r"/\*\*(.*?)\*/", text[:6000], re.S)
        if m:
            lines = [re.sub(r"^\s*\*\s?", "", ln) for ln in m.group(1).splitlines()]
            return "\n".join(lines).strip()
        lines = []
        for ln in text[:6000].splitlines():
            if ln.strip().startswith("///"):
                lines.append(ln.strip()[3:].strip())
            elif lines:
                break
        return "\n".join(lines).strip()

    def warn(self, msg, line_no):
        self.warnings.append(f"{self.rel}:{line_no + 1}: {msg}")

    # ---- the namespace-scope cursor ----
    def parse(self):
        doc = []
        while self.i < len(self.lines):
            masked, raw, dcomment = self.lines[self.i]
            s = masked.strip()
            if not s:
                if dcomment:
                    doc.append(self._doc_text(dcomment))
                else:
                    doc = []
                self.i += 1
                continue
            if s.startswith("#"):
                doc = []
                self.i += 1
                continue
            m = re.match(r"namespace\s+([\w:]+)\s*\{", s)
            if m:
                self.ns.extend(m.group(1).split("::"))
                self.i += 1
                doc = []
                continue
            if s.startswith("}"):
                if self.ns:
                    self.ns.pop()
                else:
                    self.warn("unbalanced '}' at namespace scope", self.i)
                self.i += 1
                doc = []
                continue
            template = None
            if s.startswith("template"):
                template, ok = self._take_template()
                if not ok:
                    continue
                masked, raw, dcomment = self.lines[self.i]
                s = masked.strip()
            if re.match(r"(typedef\s+)?(struct|class)\b", s):
                self._take_class(doc, template)
                doc = []
                continue
            if re.match(r"enum\b", s):
                self._take_enum(doc)
                doc = []
                continue
            if re.match(r"using\s+" + ID_RE + r"\s*=", s) or re.match(r"using\s+namespace\b", s):
                self._take_alias(doc, template)
                doc = []
                continue
            if s.startswith("static_assert") or s.startswith("extern"):
                self._skip_statement()
                doc = []
                continue
            ent = self._take_decl(doc, template, scope="ns")
            if ent:
                self.entities.append(ent)
            doc = []

    @staticmethod
    def _doc_text(comment):
        return "\n".join(re.sub(r"^\s*///<?\s?", "", ln)
                         for ln in comment.split("\n")).rstrip()

    def _take_template(self):
        """Consume a (possibly multi-line) template<...> prefix. If a
        declaration follows the closing '>' on the same line, the current line
        is replaced by that remainder (the cursor stays on it)."""
        start = self.i
        depth = 0
        seen_open = False
        parts = []
        while self.i < len(self.lines):
            masked, raw, dcomment = self.lines[self.i]
            end = None
            for j, ch in enumerate(masked):
                if ch == "<":
                    depth += 1
                    seen_open = True
                elif ch == ">":
                    depth -= 1
                    if seen_open and depth == 0:
                        end = j + 1
                        break
            if end is not None:
                parts.append(raw[:end].strip())
                if masked[end:].strip():
                    self.lines[self.i] = (masked[end:], raw[end:], dcomment)
                else:
                    self.i += 1
                return " ".join(parts), True
            parts.append(raw.strip())
            self.i += 1
        self.warn("unterminated template<>", start)
        return None, False

    def _skip_statement(self):
        while self.i < len(self.lines):
            masked, _, _ = self.lines[self.i]
            self.i += 1
            if ";" in masked:
                return

    # ---- entity parsers ----
    def _take_enum(self, doc):
        start = self.i
        masked, raw, _ = self.lines[self.i]
        m = re.match(r"\s*enum(\s+class)?\s+(" + ID_RE + ")", masked)
        if not m:
            self.warn("unparsed enum", start)
            self.i += 1
            return
        name = m.group(2)
        is_class = bool(m.group(1))
        enumerators = []
        pend = []
        # Advance to the line holding '{'; enumerators may share that line
        # (one-line enums) or follow on their own lines until '}'.
        while self.i < len(self.lines) and "{" not in self.lines[self.i][0]:
            self.i += 1
        if self.i >= len(self.lines):
            self.warn("enum without body", start)
            return
        masked, raw, dcomment = self.lines[self.i]
        seg = masked[masked.index("{") + 1:]
        while True:
            closed = "}" in seg
            if closed:
                seg = seg[:seg.index("}")]
            for tok in seg.split(","):
                tm = re.match(r"\s*(" + ID_RE + r")\s*(=.*)?$", tok.strip())
                if tm:
                    d = self._doc_text(dcomment) if dcomment else " ".join(pend)
                    enumerators.append((tm.group(1), d))
                    pend = []
            self.i += 1
            if closed:
                break
            if self.i >= len(self.lines):
                self.warn("unterminated enum", start)
                break
            masked, raw, dcomment = self.lines[self.i]
            if not masked.strip() and dcomment:
                pend.append(self._doc_text(dcomment))
                seg = ""
                continue
            seg = masked
        self.entities.append({
            "kind": "enum class" if is_class else "enum",
            "name": name, "ns": "::".join(self.ns),
            "sig": ("enum class " if is_class else "enum ") + name,
            "doc": "\n".join(doc), "line": start + 1,
            "enumerators": enumerators, "members": [],
            "template": None, "has_body": True,
        })

    def _take_class(self, doc, template, nested=False):
        """Parse struct/class. Returns the entity (also appended unless nested)."""
        start = self.i
        head_masked = []
        head_raw = []
        typedef = False
        while self.i < len(self.lines):
            masked, raw, dcomment = self.lines[self.i]
            if "{" in masked:
                k = masked.index("{")
                head_masked.append(masked[:k].strip())
                head_raw.append(raw[:k].strip())
                rest_m, rest_r = masked[k + 1:], raw[k + 1:]
                if rest_m.strip():
                    # One-line body (or members share the head line): push the
                    # remainder back as the current line for the member loop.
                    self.lines[self.i] = (rest_m, rest_r, dcomment)
                else:
                    self.i += 1
                break
            head_masked.append(masked.strip())
            head_raw.append(raw.strip())
            joined = " ".join(head_masked)
            if ";" in masked and "{" not in joined:
                self.i += 1
                return None  # forward declaration — not an entity
            self.i += 1
        joined = " ".join(head_raw)
        m = re.match(r"(typedef\s+)?(struct|class)\s+(" + ID_RE + r")?", joined)
        if not m:
            self.warn("unparsed struct/class head", start)
            return None
        typedef = bool(m.group(1))
        kind = m.group(2)
        name = m.group(3) or ""
        access_public = (kind == "struct")
        members = []
        n_private = 0
        pend_doc = []
        pend_template = None
        depth = 1
        while self.i < len(self.lines) and depth > 0:
            masked, raw, dcomment = self.lines[self.i]
            s = masked.strip()
            if not s:
                if dcomment:
                    pend_doc.append(self._doc_text(dcomment))
                else:
                    pend_doc = []
                self.i += 1
                continue
            if s.startswith("}"):
                depth -= 1
                tail = s[1:].strip()
                if typedef:
                    tm = re.match(r"(" + ID_RE + r")\s*;", tail)
                    if tm:
                        name = tm.group(1) + " (typedef struct " + (name or "?") + ")"
                self.i += 1
                break
            if s.startswith("public:"):
                access_public = True
                self.i += 1
                continue
            if s.startswith("private:") or s.startswith("protected:"):
                access_public = False
                self.i += 1
                continue
            if s.startswith("template"):
                pend_template, ok = self._take_template()
                if not ok:
                    continue
                continue
            if re.match(r"(typedef\s+)?(struct|class)\b", s):
                sub = self._take_class(pend_doc, pend_template, nested=True)
                if sub is not None:
                    if access_public:
                        members.append(sub)
                    else:
                        n_private += 1
                pend_doc = []
                pend_template = None
                continue
            if re.match(r"enum\b", s):
                mark = len(self.entities)
                self._take_enum(pend_doc)
                if len(self.entities) > mark:
                    sub = self.entities.pop()
                    if access_public:
                        members.append(sub)
                    else:
                        n_private += 1
                pend_doc = []
                continue
            if re.match(r"using\s+" + ID_RE + r"\s*=", s):
                mark = len(self.entities)
                self._take_alias(pend_doc, pend_template)
                if len(self.entities) > mark:
                    sub = self.entities.pop()
                    if access_public:
                        members.append(sub)
                    else:
                        n_private += 1
                pend_doc = []
                pend_template = None
                continue
            ent = self._take_decl(pend_doc, pend_template, scope="class")
            if ent:
                if access_public:
                    members.append(ent)
                else:
                    n_private += 1
            pend_doc = []
            pend_template = None
        ent = {
            "kind": kind, "name": name, "ns": "::".join(self.ns),
            "sig": (template + "\n" if template else "") + kind + " " + name,
            "doc": "\n".join(doc), "line": start + 1,
            "members": members, "n_private": n_private,
            "enumerators": [], "template": template, "has_body": True,
        }
        if not nested:
            self.entities.append(ent)
        return ent

    def _take_alias(self, doc, template):
        start = self.i
        text = []
        while self.i < len(self.lines):
            masked, raw, _ = self.lines[self.i]
            text.append(raw.strip())
            self.i += 1
            if ";" in masked:
                break
        joined = " ".join(text)
        m = re.match(r"using\s+(" + ID_RE + ")", joined)
        if not m:
            return  # using namespace — consumed, not an entity
        self.entities.append({
            "kind": "alias", "name": m.group(1), "ns": "::".join(self.ns),
            "sig": (template + "\n" if template else "") + joined.rstrip(";"),
            "doc": "\n".join(doc), "line": start + 1,
            "enumerators": [], "members": [], "template": template,
            "has_body": False,
        })

    def _take_decl(self, doc, template, scope):
        """Parse a function / variable / data-member declaration."""
        start = self.i
        parts_masked = []
        parts_raw = []
        paren = angle = brace = 0
        saw_eq_top = False
        first_paren_before_eq = False
        terminator = None
        guard = 0
        trailing_doc = None
        # The guard only bounds runaway misparses; generated tables (the 678-row
        # nutation initializer) legitimately span hundreds of lines.
        while self.i < len(self.lines) and guard < 5000:
            masked, raw, dcomment = self.lines[self.i]
            if trailing_doc is None and dcomment and dcomment.startswith("///<"):
                trailing_doc = self._doc_text(dcomment)
            j = 0
            cut = None
            while j < len(masked):
                ch = masked[j]
                nxt = masked[j + 1] if j + 1 < len(masked) else ""
                if ch == "o" and masked[j:j + 8] == "operator" and \
                        (j == 0 or not (masked[j - 1].isalnum() or masked[j - 1] == "_")):
                    j += 8
                    while j < len(masked) and masked[j] in "+-*/%^&|~!<>=":
                        j += 1
                    if masked[j:j + 2] in ("()", "[]"):
                        j += 2
                    continue
                if ch == "-" and nxt == ">":
                    j += 2
                    continue
                if ch == "(":
                    if paren == 0 and angle == 0 and brace == 0 and not saw_eq_top \
                            and not first_paren_before_eq:
                        first_paren_before_eq = True
                    paren += 1
                elif ch == ")":
                    paren -= 1
                elif ch == "<":
                    if nxt in "<=":
                        j += 2
                        continue
                    angle += 1
                elif ch == ">":
                    if nxt == "=":
                        j += 2
                        continue
                    if angle > 0:
                        angle -= 1
                elif ch == "=" and paren == 0 and angle == 0 and brace == 0:
                    if nxt != "=" and (j == 0 or masked[j - 1] not in "<>!+-*/%&|^="):
                        saw_eq_top = True
                elif ch == "{":
                    if paren == 0 and angle == 0 and brace == 0 and not saw_eq_top:
                        terminator = "{"
                        cut = j
                        break
                    brace += 1
                elif ch == "}":
                    brace -= 1
                elif ch == ";" and paren == 0 and brace == 0:
                    terminator = ";"
                    cut = j
                    break
                j += 1
            if terminator:
                parts_masked.append(masked[:cut])
                parts_raw.append(raw[:cut])
                if terminator == "{":
                    self._skip_body(masked, cut)
                else:
                    rest_m, rest_r = masked[cut + 1:], raw[cut + 1:]
                    if rest_m.strip():
                        # More content after ';' on this line (one-line POD
                        # members, a closing '};'): push it back.
                        self.lines[self.i] = (rest_m, rest_r, dcomment)
                    else:
                        self.i += 1
                break
            parts_masked.append(masked)
            parts_raw.append(raw)
            self.i += 1
            guard += 1
        if not terminator:
            self.warn("declaration without terminator", start)
            self.i = start + 1
            return None
        sig_raw = " ".join(p.strip() for p in parts_raw if p.strip())
        sig_masked = " ".join(p.strip() for p in parts_masked if p.strip())
        if not sig_raw:
            return None
        is_fn = first_paren_before_eq
        if is_fn:
            name = self._fn_name(sig_masked)
            sig_raw = self._trim_ctor_init(sig_raw, sig_masked)
            kind = "function"
        else:
            name = self._var_name(sig_masked)
            kind = "constant" if scope == "ns" else "data"
            sig_raw = re.sub(r"\s*=.*$", "", sig_raw).strip() or sig_raw
        if not name:
            self.warn(f"could not name declaration: {sig_raw[:80]!r}", start)
            return None
        d = "\n".join(doc) if doc else (trailing_doc or "")
        if not d and not is_fn:
            # Fall back to a trailing plain // comment on the first line.
            cm = self.lines[start][2]
            if cm:
                d = self._doc_text(cm)
        return {
            "kind": kind, "name": name, "ns": "::".join(self.ns),
            "sig": (template + "\n" if template else "") + sig_raw,
            "doc": d, "line": start + 1,
            "enumerators": [], "members": [], "template": template,
            "has_body": terminator == "{",
        }

    def _skip_body(self, masked_first, brace_pos):
        """Skip a balanced { ... } region; content after the closing brace on
        its line is pushed back as the current line."""
        depth = 0
        masked = masked_first
        pos = brace_pos
        while True:
            close = None
            for k in range(pos, len(masked)):
                ch = masked[k]
                if ch == "{":
                    depth += 1
                elif ch == "}":
                    depth -= 1
                    if depth == 0:
                        close = k
                        break
            if close is not None:
                rest_m = masked[close + 1:]
                if rest_m.strip():
                    raw = self.lines[self.i][1]
                    self.lines[self.i] = (rest_m, raw[close + 1:],
                                          self.lines[self.i][2])
                else:
                    self.i += 1
                return
            self.i += 1
            if self.i >= len(self.lines):
                return
            masked = self.lines[self.i][0]
            pos = 0

    @staticmethod
    def _fn_name(sig_masked):
        """The identifier (or operator token) before the first top-level '('."""
        paren = angle = 0
        j = 0
        last_id_end = -1
        n = len(sig_masked)
        while j < n:
            ch = sig_masked[j]
            nxt = sig_masked[j + 1] if j + 1 < n else ""
            if sig_masked[j:j + 8] == "operator" and \
                    (j == 0 or not (sig_masked[j - 1].isalnum() or sig_masked[j - 1] == "_")):
                k = j + 8
                while k < n and sig_masked[k] in "+-*/%^&|~!<>=":
                    k += 1
                if sig_masked[k:k + 2] in ("()", "[]"):
                    k += 2
                if paren == 0 and angle == 0:
                    return sig_masked[j:k].replace(" ", "")
                j = k
                continue
            if ch == "(":
                if paren == 0 and angle == 0:
                    m = re.search(r"(~?" + ID_RE + r")\s*$", sig_masked[:j])
                    return m.group(1) if m else None
                paren += 1
            elif ch == ")":
                paren -= 1
            elif ch == "<":
                if nxt in "<=":
                    j += 2
                    continue
                angle += 1
            elif ch == ">":
                if angle > 0:
                    angle -= 1
            elif ch == "-" and nxt == ">":
                j += 2
                continue
            j += 1
        return None

    @staticmethod
    def _var_name(sig_masked):
        """Declarator name(s); C-style comma declarators yield a joined list."""
        s = re.sub(r"=.*$", "", sig_masked)
        s = re.sub(r"\[[^\]]*\]", "", s)
        if "," in s and "<" not in s and "(" not in s:
            names = []
            for part in s.split(","):
                ids = re.findall(ID_RE, part)
                if ids:
                    names.append(ids[-1])
            if len(names) > 1:
                return ", ".join(names)
        m = re.findall(ID_RE, s)
        return m[-1] if m else None

    @staticmethod
    def _trim_ctor_init(sig_raw, sig_masked):
        """Cut a constructor's member-init list (': x_(...), ...') if present."""
        paren = angle = 0
        depth_hit_zero_after_open = False
        opened = False
        j = 0
        n = len(sig_masked)
        close = None
        while j < n:
            ch = sig_masked[j]
            nxt = sig_masked[j + 1] if j + 1 < n else ""
            if ch == "(":
                if paren == 0 and angle == 0:
                    opened = True
                paren += 1
            elif ch == ")":
                paren -= 1
                if opened and paren == 0 and angle == 0:
                    close = j
                    break
            elif ch == "<":
                if nxt in "<=":
                    j += 2
                    continue
                angle += 1
            elif ch == ">":
                if angle > 0:
                    angle -= 1
            elif ch == "-" and nxt == ">":
                j += 2
                continue
            j += 1
        if close is None:
            return sig_raw
        tail = sig_masked[close + 1:]
        m = re.search(r"(?<!:):(?!:)", tail)
        if m:
            keep = close + 1 + m.start()
            return sig_raw[:keep].rstrip()
        return sig_raw


def _inherit_overload_docs(ents):
    """A bare overload immediately following a documented same-name function
    shares that doc (the C++ overload-set convention)."""
    prev = None
    for e in ents:
        if e["kind"] == "function" and not e["doc"].strip() and prev is not None \
                and prev["kind"] == "function" and prev["name"] == e["name"] \
                and prev["ns"] == e["ns"] and prev["doc"].strip():
            e["doc"] = prev["doc"]
        prev = e
        if e["members"]:
            _inherit_overload_docs(e["members"])


def parse_header(path, rel):
    p = HeaderParser(path, rel)
    p.parse()
    _inherit_overload_docs(p.entities)
    return p


# ============================ doc-block helpers ==============================

def doc_block(path):
    """Extract the leading /** ... */ or /// doc-block of a header."""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read(6000)
    except OSError:
        return ""
    m = re.search(r"/\*\*(.*?)\*/", text, re.S)
    if m:
        lines = [re.sub(r"^\s*\*\s?", "", ln) for ln in m.group(1).splitlines()]
    else:
        lines = []
        for ln in text.splitlines():
            if ln.strip().startswith("///"):
                lines.append(ln.strip()[3:].strip())
            elif lines:
                break
    return "\n".join(lines).strip()


def linkify(text):
    """Escape, then link theory notes and shorten file refs."""
    esc = html.escape(text)
    esc = re.sub(r"design/derivations/([\w\-\.]+\.md)",
                 r'<a href="' + DESIGN_BASE + r'\1">design/derivations/\1</a>', esc)
    return esc


def gates_for(text):
    return sorted(set(re.findall(
        r"\b(OR1|TB1|SRP1|ATM1|FM1|RK1|NUT1|EX1|EX2|DOC1|GAL1|GEOPOT|FRAME1|FRAME2|"
        r"TIME1|SC1|EPH|DRAG1|AD1|DS1|D1|D2|CR1B?|C1|C2|E1|E2|E3|F1|F2|F3|G1|H1|INJ1|BUG1|"
        r"W1[1-9]|W2[0-4]|B[1-3])\b", text)))


def first_para(text, limit=1200):
    p = text.split("\n\n")[0] if text else ""
    return p[:limit] + (" …" if len(p) > limit else "")


def first_sentence(text, limit=240):
    p = text.replace("\n", " ").strip()
    m = re.match(r"(.+?\.)\s", p + " ")
    s = m.group(1) if m else p
    return s[:limit] + (" …" if len(s) > limit else "")


def anchor_for(ent):
    base = re.sub(r"\W+", "-", ent["name"]).strip("-")
    return f"{base}-L{ent['line']}"


# ============================ page renderers =================================

def render_entity(ent, mod):
    a = anchor_for(ent)
    kind = ent["kind"]
    sig = html.escape(ent["sig"])
    body = [f'<div class="ent" id="{a}">'
            f'<span class="kind">{html.escape(kind)}</span>'
            f'<code>{html.escape((ent["ns"] + "::" if ent["ns"] else "") + ent["name"])}</code>'
            f'<pre class="sig">{sig}</pre>']
    if ent["doc"]:
        body.append(f'<div class="entdoc">{linkify(first_para(ent["doc"]))}</div>')
    if ent["enumerators"]:
        rows = "".join(f"<tr><td><code>{html.escape(n)}</code></td>"
                       f"<td>{linkify(d)}</td></tr>" for n, d in ent["enumerators"])
        body.append('<table class="member-tbl"><tr><th>Enumerator</th><th>Meaning</th></tr>'
                    + rows + "</table>")
    if ent["members"]:
        rows = []
        for mem in ent["members"]:
            msig = html.escape(mem["sig"])
            rows.append(f"<tr><td><code>{html.escape(mem['name'])}</code><br>"
                        f"<span class='muted'>{html.escape(mem['kind'])}</span></td>"
                        f"<td><pre class='sig' style='margin:0'>{msig}</pre>"
                        f"<div class='entdoc'>{linkify(first_sentence(mem['doc']))}</div></td></tr>")
        body.append('<table class="member-tbl"><tr><th>Member</th><th>Signature</th></tr>'
                    + "".join(rows) + "</table>")
    if ent.get("n_private"):
        body.append(f'<p class="muted">+ {ent["n_private"]} private member(s) — '
                    'implementation detail, not part of the public surface.</p>')
    body.append("</div>")
    return "".join(body)


def dedupe_prototypes(ents):
    """Drop bodiless function declarations shadowed by a same-name definition
    (the verbatim-port prototype block)."""
    defined = {(e["ns"], e["name"]) for e in ents
               if e["kind"] == "function" and e["has_body"]}
    out = []
    for e in ents:
        if e["kind"] == "function" and not e["has_body"] \
                and (e["ns"], e["name"]) in defined:
            continue
        out.append(e)
    return out


def parse_gate_registry():
    """The registered gates of run_acceptance.ps1: {id: (target, kind,
    registered_comment)}. Includes the dynamic I0 foundation gates (the
    $modules sweep and the special I0.sgp4 verification), mirrored
    mechanically from the same script."""
    registry = {}
    comment = []
    text = ""
    with open(os.path.join(ROOT, "tools", "run_acceptance.ps1"), encoding="utf-8",
              errors="replace") as f:
        for ln in f:
            text += ln
            s = ln.strip()
            if s.startswith("#"):
                comment.append(s.lstrip("# "))
                continue
            m = re.match(r"(ExeGate|ScriptGate|VerifierGate|GrepGate)\s+'([^']+)'\s+'([^']+)'", s)
            if m:
                kind, gid, name = m.groups()
                registry[gid] = (name, kind, " ".join(comment)[:400])
            comment = []
    m = re.search(r"\$modules\s*=\s*((?:'[\w]+'[,\s]*)+)", text)
    if m:
        for name in re.findall(r"'([\w]+)'", m.group(1)):
            registry[f"I0.{name}"] = (name, "FoundationGate",
                                      "foundation sweep: the module test exits 0 "
                                      "on every validation")
    registry["I0.sgp4"] = ("test_sgp4", "FoundationGate",
                           "the frozen invariant: the official SGP4-VER suite — "
                           "33/33 satellites, 623/623 points, asserted verbatim "
                           "from the test output on every validation")
    return registry


def theory_note_meta(fn):
    """(title, lead) of a design/derivations note."""
    try:
        with open(os.path.join(DERIV, fn), encoding="utf-8", errors="replace") as f:
            text = f.read(3000)
    except OSError:
        return "", ""
    title = ""
    lead = ""
    for ln in text.splitlines():
        if ln.startswith("# ") and not title:
            title = ln[2:].strip()
        elif title and ln.strip() and not ln.startswith("#"):
            lead = ln.strip()
            break
    return title, lead


def notes_for_module(mod):
    """Theory notes referenced anywhere in the module's headers (mechanical)."""
    moddir = os.path.join(SRC, mod)
    notes = set()
    for fn in sorted(os.listdir(moddir)):
        if not fn.endswith(".h") and not fn.endswith(".cpp"):
            continue
        try:
            with open(os.path.join(moddir, fn), encoding="utf-8", errors="replace") as f:
                text = f.read()
        except OSError:
            continue
        notes.update(re.findall(r"design/derivations/([\w\-\.]+\.md)", text))
    return sorted(n for n in notes if os.path.exists(os.path.join(DERIV, n)))


def gate_module_map(registry):
    """ExeGate id → the src/ modules its test's main.cpp DIRECTLY includes
    (mechanical: a gate exercises the modules it compiles against)."""
    out = {}
    for gid, (name, kind, _) in registry.items():
        if kind not in ("ExeGate", "FoundationGate"):
            continue
        main = os.path.join(ROOT, "tests", name, "main.cpp")
        if not os.path.exists(main):
            continue  # the example gates live in examples/
        try:
            with open(main, encoding="utf-8", errors="replace") as f:
                text = f.read()
        except OSError:
            continue
        mods = set(re.findall(r'#include\s+"(?:\.\./)*(?:src/)?([a-z_0-9]+)/[\w\.]+\.h"',
                              text))
        out[gid] = {m for m in mods if os.path.isdir(os.path.join(SRC, m))}
    return out


def gates_for_module(mod, gate_mods):
    """Gates covering the module: named in its header doc-blocks, or whose
    test directly includes one of its headers (both mechanical)."""
    moddir = os.path.join(SRC, mod)
    ids = set()
    for fn in sorted(os.listdir(moddir)):
        if not fn.endswith(".h"):
            continue
        ids.update(gates_for(doc_block(os.path.join(moddir, fn))))
    ids.update(gid for gid, mods in gate_mods.items() if mod in mods)
    return sorted(ids)


def module_page(mod, blurb, parsed, registry, gate_mods, diag):
    moddir = os.path.join(SRC, mod)

    diagrams_sec = "".join(diag[d] for d in MODULE_DIAGRAMS.get(mod, []))

    tier_sec = ""
    if mod in MODULE_TIER_NOTE:
        tier_sec = (f'<div class="card"><p><b>{html.escape(MODULE_TIER_NOTE[mod])}</b>'
                    "</p></div>")

    # --- Theory ---
    notes = notes_for_module(mod)
    theory_rows = []
    for n in notes:
        title, lead = theory_note_meta(n)
        theory_rows.append(f"<tr><td><a href='{DESIGN_BASE}{n}'>{n}</a></td>"
                           f"<td>{html.escape(title)}</td>"
                           f"<td>{html.escape(lead[:240])}</td></tr>")
    theory = ("<h2>Theory</h2>"
              "<p>The notes this module's headers cite — the derivations the code "
              "fell out of (theory-first; see <a href='theory.html'>the full index</a>).</p>"
              "<table><tr><th>Note</th><th>Title</th><th>Lead</th></tr>"
              + "".join(theory_rows) + "</table>") if theory_rows else (
              "<h2>Theory</h2><p class='muted'>No design/derivations note is cited by "
              "this module's headers; the doc-blocks below carry the basis.</p>")

    # --- Usage ---
    usage_parts = ["<h2>Usage</h2>"]
    if mod in MODULE_USAGE_NOTE:
        usage_parts.append(f"<p>{html.escape(MODULE_USAGE_NOTE[mod])}</p>")
    for src_file, tag, gate in MODULE_USAGE.get(mod, []):
        code = html.escape(extract_snippet(os.path.join(ROOT, src_file), tag))
        usage_parts.append(
            f"<pre>{code}</pre>"
            f'<p class="provenance">Extracted verbatim from <code>{src_file}</code> — '
            f"compiled and run on every validation by gate <b>{gate}</b>.</p>")
    if not MODULE_USAGE.get(mod) and mod not in MODULE_USAGE_NOTE:
        usage_parts.append("<p class='muted'>No standalone usage snippet; see the "
                           "<a href='guide.html'>guide</a>.</p>")
    usage = "".join(usage_parts)

    # --- Gates & measured grades ---
    gids = gates_for_module(mod, gate_mods)
    gate_rows = []
    for gid in gids:
        if gid in registry:
            name, kind, note = registry[gid]
            gate_rows.append(f"<tr><td><b>{gid}</b></td><td><code>{html.escape(name)}</code></td>"
                             f"<td>{html.escape(note)}</td></tr>")
        else:
            gate_rows.append(f"<tr><td><b>{gid}</b></td><td colspan='2' class='muted'>"
                             "referenced in doc-blocks; not a statically registered "
                             "gate line (a check id or a dynamic gate)</td></tr>")
    tests_sec = ("<h2>Test cases — gates and measured grades</h2>"
                 "<p>The acceptance gates covering this module — named in its header "
                 "doc-blocks, or whose test directly includes one of its headers. The "
                 "registered note beside each gate records what it asserts and the "
                 "MEASURED grade (the no-perceived-fidelity rule: stated grades are "
                 "measured, never aspirational). The <a href='tests.html'>full "
                 "registry</a> lists every gate.</p>"
                 "<table><tr><th>Gate</th><th>Target</th><th>What it asserts (registered "
                 "note, measured grades)</th></tr>" + "".join(gate_rows) + "</table>"
                 ) if gate_rows else (
                 "<h2>Test cases</h2><p class='muted'>No gate ids in this module's "
                 "doc-blocks; coverage arrives via its consumers' gates.</p>")

    # --- References ---
    refs = MODULE_REFS.get(mod, [])
    refs_sec = ("<h2>References</h2><ul>"
                + "".join(f"<li>{linkify(r)}</li>" for r in refs)
                + "</ul>") if refs else ""

    # --- Per-header API reference (Q1) ---
    sections = []
    toc = []
    for fn in sorted(os.listdir(moddir)):
        if not fn.endswith(".h") and not fn.endswith(".cpp"):
            continue
        block = doc_block(os.path.join(moddir, fn))
        g = ", ".join(gates_for(block)) or "—"
        hid = "hdr-" + re.sub(r"\W+", "-", fn)
        toc.append(f'<a href="#{hid}"><code>{fn}</code></a>')
        sec = [f'<h3 id="{hid}"><code>src/{mod}/{fn}</code></h3>',
               f'<p class="muted">Gates referenced in the doc-block: {g}</p>',
               f'<div class="card"><pre class="docblock">{linkify(first_para(block, 2200))}</pre></div>']
        p = parsed.get((mod, fn))
        if p:
            pub = [e for e in dedupe_prototypes(p.entities)
                   if not (set(e["ns"].split("::")) & INTERNAL_NAMESPACES)]
            intern = [e for e in dedupe_prototypes(p.entities)
                      if set(e["ns"].split("::")) & INTERNAL_NAMESPACES]
            if pub:
                sec.append("".join(render_entity(e, mod) for e in pub))
            if intern:
                inner = "".join(render_entity(e, mod) for e in intern)
                sec.append(f"<details><summary>internal namespace — {len(intern)} "
                           "entities (implementation detail; not part of the public "
                           "surface)</summary>" + inner + "</details>")
            if not pub and not intern:
                sec.append('<p class="muted">No namespace-scope declarations '
                           '(an umbrella or table include).</p>')
        sections.append("".join(sec))
    api_sec = ("<h2>API reference</h2>"
               '<p class="muted">Headers: ' + " · ".join(toc) + "</p>"
               + "".join(sections))

    body = (f"<p>{html.escape(blurb)}</p>" + tier_sec + diagrams_sec
            + theory + usage + tests_sec + refs_sec + api_sec)
    return page(f"src/{mod}/ — module documentation", body)


def api_index_page(parsed):
    rows = []
    n_fn = n_ty = 0
    flat = []
    for (mod, fn), p in parsed.items():
        for e in dedupe_prototypes(p.entities):
            if set(e["ns"].split("::")) & INTERNAL_NAMESPACES:
                continue
            flat.append(((e["name"].lower(), e["ns"], mod, fn, e["line"]), e, mod, fn))
            if e["kind"] == "function":
                n_fn += 1
            else:
                n_ty += 1
            n_fn += sum(1 for m2 in e["members"] if m2["kind"] == "function")
    flat.sort(key=lambda t: t[0])
    for _, e, mod, fn in flat:
        a = anchor_for(e)
        extra = ""
        if e["members"]:
            pubfn = sum(1 for m in e["members"] if m["kind"] == "function")
            if pubfn:
                extra = f' <span class="muted">({pubfn} public member function(s))</span>'
        rows.append(
            f"<tr><td><a href='module_{mod}.html#{a}'><code>{html.escape(e['name'])}</code></a>{extra}</td>"
            f"<td>{html.escape(e['kind'])}</td>"
            f"<td><code>{html.escape(e['ns'])}</code></td>"
            f"<td><code>src/{mod}/{fn}</code></td>"
            f"<td>{linkify(first_sentence(e['doc']))}</td></tr>")
    body = (
        "<p>Every namespace-scope entity in the public headers, parsed from the source "
        "by <code>tools/gen_docs.py</code> (strict scanner — an unclassifiable construct "
        "fails generation, so this index is complete by construction). Click a name for "
        "its full signature and, for classes, every public member signature. Internal "
        "namespaces (<code>detail</code>, <code>msis</code>) are enumerated on the module "
        "pages under a collapsed section.</p>"
        f"<p class='muted'>{len(rows)} public entities ({n_fn} functions including public "
        f"members, {n_ty} types/aliases/constants).</p>"
        "<table><tr><th>Name</th><th>Kind</th><th>Namespace</th><th>Header</th>"
        "<th>Summary</th></tr>" + "".join(rows) + "</table>")
    return page("API index", body)


# ============================ the help guide =================================

def extract_snippet(path, tag):
    """Extract the region between '// [guide:tag]' and '// [guide:end]'."""
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    out = []
    grab = False
    for ln in lines:
        if re.search(r"//\s*\[guide:" + re.escape(tag) + r"\]", ln):
            grab = True
            continue
        if grab and re.search(r"//\s*\[guide:end\]", ln):
            # Dedent by the common leading whitespace.
            keep = [l for l in out if l.strip()]
            if keep:
                ind = min(len(l) - len(l.lstrip()) for l in keep)
                out = [l[ind:] if len(l) > ind else l for l in out]
            return "\n".join(out).rstrip()
        if grab:
            out.append(ln)
    raise SystemExit(f"gen_docs: missing guide marker [guide:{tag}] in {path}")


def guide_page(parsed, diag):
    """The task-oriented help guide. Snippets are extracted verbatim from the
    gate-compiled examples; the symbols each section relies on are verified to
    exist in the parsed API model (generation fails otherwise)."""
    model_names = set()
    for p in parsed.values():
        for e in p.entities:
            model_names.add(e["name"])
            for mem in e["members"]:
                model_names.add(mem["name"])

    def require(*names):
        missing = [n for n in names if n not in model_names]
        if missing:
            raise SystemExit(f"gen_docs: guide references unknown API symbols: {missing}")

    qs = os.path.join(EXAMPLES, "quickstart.cpp")
    gt = os.path.join(EXAMPLES, "ground_track.cpp")

    def snip(path, tag, gate):
        code = html.escape(extract_snippet(path, tag))
        rel = os.path.relpath(path, ROOT).replace("\\", "/")
        return (f"<pre>{code}</pre>"
                f'<p class="provenance">Extracted verbatim from <code>{rel}</code> — '
                f"compiled and run on every validation by gate <b>{gate}</b>.</p>")

    require("Propagator", "DqSgp4Propagator", "StateVector", "State", "TleData",
            "TleElements", "ModelSelector", "ModelConfiguration", "propagate",
            "from_tle", "DqForceOptions", "PropagatorMode", "propagate_adaptive",
            "AdaptiveResult", "TrackedValue", "ThreeErrors", "gcrs_to_itrs",
            "make_third_body_force", "make_srp_force", "make_drag")

    body = f"""
<p>Task-oriented recipes for the library's main jobs. Every snippet below is extracted
verbatim from a source file that the acceptance suite compiles and runs (gates EX1/EX2),
so the guide cannot drift from working code. The <a href="api.html">API index</a> has
every signature; <a href="index.html">the overview</a> tabulates what each component is verified against.</p>

<h2 id="both-ways">1. Propagate a TLE with either propagator</h2>
{diag["sgp4_vs_dqsgp4"]}
<p>The numerical propagator (<code>dynamics::DqSgp4Propagator</code>) integrates
dual-quaternion dynamics; the analytical SGP4/SDP4
(<code>sgp4::Propagator</code>, verified against the published test set of 33
satellites and 623 points) implements the model that defines a TLE's orbital
elements. Because of that, the numerical propagator recovers its initial state by
evaluating the analytical model at the epoch — short-period corrections at the
kilometre scale would be lost by reading the elements as osculating values. The
analytical headers are not part of the umbrella, so the snippet includes them
explicitly. Parse once, then:</p>
{snip(qs, "sgp4", "EX2")}
<p>The DQ facade exposes the same <code>propagate(minutes)</code> verb
(<code>dynamics::Propagatable</code> is the shared concept; gate F2 runs one generic
driver over BOTH families and asserts bit-equality with direct calls):</p>
{snip(qs, "facade", "EX2")}

<h2 id="forces">2. Add perturbation forces — presets and extra_forces</h2>
<p>The default DQ force model is gravitational only (the unified geopotential —
monopole + zonal in one Cunningham pass, gate GEOPOT). Perturbations are opt-in two
ways. The <code>DqForceOptions</code> presets assemble independently verified forces from
the TLE epoch — lunisolar third-body (gate TB1, JPL DE430-verified ephemeris),
Vallado 8-4 table drag (gate ATM1), cannonball SRP (gate SRP1). Gate FM1 measures the
presets on SGP4-VER satellite 00005 over 60 minutes: lunisolar moves the endpoint
12.6&nbsp;m, SRP 1.39&nbsp;m, drag 0.33&nbsp;m — and asserts the default stays
bit-identical to the gravitational-only model:</p>
{snip(qs, "presets", "EX2")}
<p>Arbitrary forces append after the gravitational core via the explicit constructor's
<code>extra_forces</code> (each a <code>(State, ConstantsProvider) → Wrench</code>
lambda; see <code>forces::make_drag</code>, <code>forces::make_third_body_force</code>,
<code>forces::make_srp_force</code> for ready-made builders; a space-weather
density model plugs the same <code>DensityModel</code> interface).</p>

<h2 id="adaptive">3. Adaptive stepping</h2>
<p><code>dynamics::Propagator::propagate_adaptive</code> runs the RKF7(8) embedded
pair with per-step error control instead of the fixed cadence of
<code>propagate_to</code>. Gate AD1 asserts it is bit-identical to the standalone
adaptive loop and exercises the accept/reject diagnostics it returns:</p>
{snip(qs, "adaptive", "EX2")}

<h2 id="itrs">4. Earth-fixed output — ITRS and ground tracks</h2>
{diag["frame_chain"]}
<p>SGP4 and the DQ propagator emit TEME. The Earth-fixed chain is
r<sub>ITRS</sub> = W·R3(GAST)·N·P·r<sub>GCRS</sub>:
<code>astronomy::gcrs_to_itrs</code> (IAU 2000A nutation from the in-repo SOFA table,
GMST06/GAST, polar motion). Gate NUT1 measures the chain against erfa
<code>c2t06a</code> at ~0.4&nbsp;mas (polar motion bit-exact against
<code>pom00</code>). The gated ground-track example walks the whole pipeline —
TLE → DQSGP4 → TEME → GCRS → ITRS → subpoint:</p>
{snip(gt, "itrs", "EX1")}
<p>Demo honesty (the example prints real grades): using the TLE's UTC epoch for TT
shifts the subpoint longitude ~0.29°; zero polar motion ≤ ~9&nbsp;m; the subpoint is
geocentric (geodetic conversion adds up to ~0.19° at mid-latitudes).</p>

<h2 id="errors">5. Reading the three-error budget</h2>
<p>Every computed quantity is a <code>math::TrackedValue&lt;T&gt;</code>: a value
plus a <code>ThreeErrors</code> budget — <code>measurement</code> (physical input
σ), <code>precision</code> (representation/rounding), <code>accuracy</code> (model
truncation), each propagated through every operation (gate W12 and the per-module
gates assert the propagation). Read them directly:</p>
{snip(qs, "errors", "EX2")}
<p>The channels answer different questions: precision tells you whether a wider
<code>T</code> would change the digits; accuracy tells you what the MODEL cannot
claim regardless of arithmetic; measurement carries published σ from constants
(e.g. GM☉). <code>total()</code> is the conservative triangle-inequality bound.</p>
<p>A note on the accuracy bound of a numerically propagated state. The fixed-step
integrators record a conservative truncation bound on each step, assigned by
component: a position-scale bound of (h²/2)·A on the translation-carrying part of
the pose, a velocity-scale bound of h·A on the twist, and the corresponding
rotational bounds — which are exactly zero when no torque acts — on the rotational
components (the derivation is in
<code>runge_kutta_lie_group.md</code> §6.2, and the test suite checks the
assignment). For a single 30-second RK4 step on a low-Earth orbit this works out to
about 4&nbsp;km of claimed position accuracy — deliberately conservative, since
only the acceleration magnitude is available to bound the higher derivatives. Over
many fixed steps the accumulated bound eventually exceeds the state itself, at
which point the propagated bound legitimately diverges; it cannot certify a long
arc. For a meaningful accuracy figure over a long arc, use
<code>propagate_adaptive</code>, whose embedded 7(8) difference measures the actual
per-step error, and read the model accuracy of each force from its own
budget.</p>

<h2 id="choosing-t">6. Choosing T — double vs cpp_bin_float_50</h2>
<p>Everything is templated on the numeric type. <code>double</code> is the working
default; <code>boost::multiprecision::cpp_bin_float_50</code> tightens the PRECISION
channel (the per-module gates measure this — e.g. gate SC1 shows series precision
scaling with T while accuracy stays pinned by truncation). Wider T does NOT improve
accuracy: an empirical model keeps its published grade regardless of the
arithmetic width (more digits cannot improve a fitted model).
Instantiate with the type you need:</p>
{snip(qs, "precision", "EX2")}

<h2 id="more">Where to go next</h2>
<ul>
<li><a href="api.html">API index</a> — every public entity and signature.</li>
<li><a href="theory.html">Theory notes</a> — the derivations each module falls out of.</li>
<li><a href="tests.html">Acceptance gates</a> — what is measured, where.</li>
<li><code>README.md</code> — build instructions (MSBuild on sgp4.sln) and quick starts.</li>
</ul>
"""
    return page("Help guide", body)


# ============================ stock pages ====================================

def theory_page():
    rows = []
    for fn in sorted(os.listdir(DERIV)):
        if not fn.endswith(".md"):
            continue
        try:
            with open(os.path.join(DERIV, fn), encoding="utf-8", errors="replace") as f:
                text = f.read(3000)
        except OSError:
            continue
        title = ""
        lead = ""
        for ln in text.splitlines():
            if ln.startswith("# ") and not title:
                title = ln[2:].strip()
            elif title and ln.strip() and not ln.startswith("#"):
                lead = ln.strip()
                break
        rows.append(f"<tr><td><a href='{DESIGN_BASE}{fn}'>{fn}</a></td>"
                    f"<td>{html.escape(title)}</td><td>{html.escape(lead[:300])}</td></tr>")
    body = ("<p>The theory notes are the substance of the theory-first method: every module's "
            "exhaustive derivation precedes its code, and every fidelity claim names its "
            "independent oracle and tolerance.</p>"
            "<table><tr><th>Note</th><th>Title</th><th>Lead</th></tr>"
            + "".join(rows) + "</table>")
    return page("Theory notes (design/derivations)", body)


def tests_page(registry):
    rows = []
    for gid, (name, kind, note) in registry.items():
        rows.append(f"<tr><td><b>{gid}</b></td><td><code>{html.escape(name)}</code></td>"
                    f"<td>{kind}</td><td>{html.escape(note)}</td></tr>")
    body = ("<p>The acceptance suite (<code>tools/run_acceptance.ps1</code>): every gate runs "
            "on every validation; the fixed reference points (the bit-exact SGP4 regression, "
            "test_sgp4 33/33 vs SGP4-VER) hold at every commit. The table lists the "
            "registered gate lines plus the I0 foundation sweep (mirrored from the same "
            "script); the remaining I0 checks (build, golden files, the Octave/Python "
            "verifier sweeps) are assembled dynamically and appear in the console "
            "report.</p>"
            "<table><tr><th>Gate</th><th>Target</th><th>Kind</th><th>Registered note</th></tr>"
            + "".join(rows) + "</table>")
    return page("Acceptance gates", body)


def index_page(diag):
    oracle_rows = [
        ("Time-scale conversions / epochs", "IERS constants; round-trip identities; the in-repo SGP4 epoch handling (TIME1)"),
        ("Frames: precession, obliquity, frame products", "ERFA (IAU SOFA) bit/identity gates (FRAME1/FRAME2)"),
        ("Sun/Moon ephemerides (Meeus series)", "JPL DE430, the in-repo numerical ephemeris — NOT the same analytical theory (EPH/FRAME2)"),
        ("Nutation + GAST + polar motion + GCRS→ITRS", "erfa nut06a / gst06a / pom00 / c2t06a; polar motion bit-exact, chain ~0.4 mas (NUT1)"),
        ("Geopotential (monopole+zonal+tesseral)", "summed legacy forces + closed-form J₂ + the monopole identity (GEOPOT, D2)"),
        ("EGM2008 zonal provenance", "the in-repo 15-digit NGA coefficient file; cross-site double-entry (GAL1)"),
        ("Static atmosphere (Vallado 8-4)", "the in-repo born-digital table: node-exact + chain-verified (ATM1)"),
        ("Third-body (Sun/Moon tidal force)", "JPL DE430 positions; Battin-stable form vs naive at cpp_bin_float_50 (TB1)"),
        ("Solar radiation pressure", "P₁ᴬᵁ generated from L☉/(4π·AU²·c); TSI basis 1361 W/m²; geometry identities (SRP1)"),
        ("Integrators (euler/RK4/RKF78)", "exact tableau order-conditions; analytic Kepler convergence orders (RK1, W20–W22)"),
        ("SGP4/SDP4 (the analytical reference)", "the official SGP4-VER suite, 33/33 satellites, 623/623 points, plus a bit-exact regression baseline"),
        ("Documentation freshness", "regenerate-and-diff: docs/ must equal what gen_docs.py emits from the current tree (DOC1)"),
    ]
    otab = "".join(f"<tr><td>{html.escape(a)}</td><td>{html.escape(b)}</td></tr>"
                   for a, b in oracle_rows)
    mods = "".join(
        f"<tr><td><a href='module_{m}.html'><code>src/{m}/</code></a></td>"
        f"<td>{html.escape(b)}</td></tr>"
        for m, b in MODULES if os.path.isdir(os.path.join(SRC, m)))
    body = f"""
<div class="card"><p><b>DQSGP4</b> is a C++20 satellite orbit propagator built on
dual-quaternion rigid-body dynamics. The spacecraft state — orientation and position
together — is carried as a unit dual quaternion on SE(3) and advanced by screw-motion
integration (fixed-step RK4 or adaptive RKF7(8)), with the attitude represented by a
quaternion throughout. Every computed value is a
<code>math::TrackedValue&lt;T&gt;</code> carrying three separate error bounds —
measurement uncertainty, numerical precision, and model accuracy — and the numeric
type <code>T</code> is a template parameter: the same code runs in
<code>double</code> or in 50-digit <code>cpp_bin_float_50</code>. The library also
contains a complete implementation of the analytical SGP4/SDP4 theory, verified
against the published test set (33 satellites, 623 evaluation points); it serves as
the verification reference for the numerical propagator and recovers initial states
from two-line element sets.</p>
<p>One include provides the public surface: <code>#include "dqsgp4.h"</code>. Start
with the <a href="guide.html">guide</a>; the <a href="api.html">API index</a> lists
every public entity. Two complete example programs are built and run by the test
suite: <code>examples/quickstart.cpp</code> (both propagators, the perturbation
options, adaptive stepping, error budgets, and precision selection) and
<code>examples/ground_track.cpp</code> (a ground track computed from a TLE through
the numerical propagator and the Earth-orientation chain).</p></div>

<h2>How the library is verified</h2>
<p>Each module is developed from a written derivation (the
<a href="theory.html">theory notes</a>), implemented, and then verified against an
independent reference before any accuracy figure is stated. Every verified property is
registered as a test in the <a href="tests.html">acceptance suite</a>, which runs as a
whole on every change. Accuracy claims in this documentation are measured results from
those tests, not estimates.</p>

<h2>The two propagators</h2>
{diag["sgp4_vs_dqsgp4"]}

<h2>Modules</h2>
<table><tr><th>Module</th><th>What it provides</th></tr>{mods}</table>
{diag["layer_map"]}

<h2>Verification — what each component is checked against</h2>
<table><tr><th>Component</th><th>Verified against</th></tr>{otab}</table>

<h2>References</h2>
<ul>
<li>Hoots &amp; Roehrich (1980), <i>Spacetrack Report #3</i>; Vallado, Crawford, Hujsak &amp; Kelso (2006), <i>Revisiting Spacetrack Report #3</i> (AIAA 2006-6753) — the SGP4/SDP4 model and verification suite.</li>
<li>Vallado, <i>Fundamentals of Astrodynamics and Applications</i> — the exponential atmosphere table (8-4), AstroLib reference code (in-repo).</li>
<li>Meeus, <i>Astronomical Algorithms</i> (1998) — solar (§25) and lunar (§47) series.</li>
<li>IERS Conventions (2010), TN36 — frames, nutation, polar motion, EGM2008 zonals.</li>
<li>IAU SOFA / ERFA — the frame-chain oracles (pmat06, nut06a, gst06a, pom00, c2t06a).</li>
<li>JPL DE430 — the ephemeris oracle (in-repo Sun/Moon table).</li>
<li>Fehlberg (1968), NASA TR R-287 — the RKF7(8) tableau; Butcher (2003) — RK order theory.</li>
<li>NGA EGM2008 — the gravity-field coefficients (in-repo 15-digit file).</li>
</ul>
"""
    return page("DQSGP4 — library documentation", body)


# ============================ driver =========================================

def undocumented_public_symbols(parsed):
    """Public entities / public member functions without a doc-comment."""
    out = []
    for (mod, fn), p in parsed.items():
        for e in dedupe_prototypes(p.entities):
            if set(e["ns"].split("::")) & INTERNAL_NAMESPACES:
                continue
            if not e["doc"].strip():
                out.append(f"src/{mod}/{fn}: {e['kind']} {e['name']}")
            for m in e["members"]:
                if m["kind"] == "function" and not m["doc"].strip():
                    out.append(f"src/{mod}/{fn}: {e['name']}::{m['name']}")
    return out


def build_all():
    """Parse every header, render every page. Returns (pages, warnings, undoc)."""
    parsed = {}
    warnings = []
    for mod, _ in MODULES:
        moddir = os.path.join(SRC, mod)
        if not os.path.isdir(moddir):
            continue
        for fn in sorted(os.listdir(moddir)):
            if not fn.endswith(".h"):
                continue
            rel = f"src/{mod}/{fn}"
            p = parse_header(os.path.join(moddir, fn), rel)
            parsed[(mod, fn)] = p
            warnings.extend(p.warnings)

    registry = parse_gate_registry()
    gate_mods = gate_module_map(registry)
    diag = gen_diagrams.diagrams()
    pages = {"index.html": index_page(diag),
             "guide.html": guide_page(parsed, diag),
             "api.html": api_index_page(parsed),
             "theory.html": theory_page(),
             "tests.html": tests_page(registry)}
    for mod, blurb in MODULES:
        if not os.path.isdir(os.path.join(SRC, mod)):
            continue
        pages[f"module_{mod}.html"] = module_page(mod, blurb, parsed, registry,
                                                  gate_mods, diag)
    return pages, warnings, undocumented_public_symbols(parsed)


def main():
    check = "--check" in sys.argv
    pages, warnings, undoc = build_all()
    if warnings:
        print("gen_docs: PARSER WARNINGS (unclassified constructs):", file=sys.stderr)
        for w in warnings:
            print("  " + w, file=sys.stderr)
        sys.exit(1)
    # The doc-completeness ratchet: count may only go DOWN (Q4 lowers the
    # baseline alongside the source fixes; end state 0 pins full coverage).
    print(f"gen_docs: {len(undoc)} undocumented public symbols "
          f"(ratchet baseline {UNDOC_BASELINE})")
    if len(undoc) > UNDOC_BASELINE:
        print("gen_docs: DOC-COVERAGE REGRESSION — new undocumented public symbols:",
              file=sys.stderr)
        for u in undoc:
            print("  " + u, file=sys.stderr)
        sys.exit(1)
    if "--list-undocumented" in sys.argv:
        for u in undoc:
            print("  " + u)
    if check:
        drift = []
        for name, content in pages.items():
            path = os.path.join(DOCS, name)
            try:
                with open(path, encoding="utf-8") as f:
                    on_disk = f.read()
            except OSError:
                drift.append(name + " (missing)")
                continue
            if on_disk != content:
                drift.append(name)
        stale = [fn for fn in os.listdir(DOCS)
                 if fn.endswith(".html") and fn not in pages]
        if drift or stale:
            if drift:
                print("gen_docs --check: STALE (regenerate and commit): "
                      + ", ".join(drift), file=sys.stderr)
            if stale:
                print("gen_docs --check: ORPHANED pages in docs/: "
                      + ", ".join(stale), file=sys.stderr)
            sys.exit(1)
        print(f"gen_docs --check: docs/ is fresh ({len(pages)} pages).")
        return
    os.makedirs(DOCS, exist_ok=True)
    for name, content in pages.items():
        with open(os.path.join(DOCS, name), "w", encoding="utf-8") as f:
            f.write(content)
    print("wrote docs/: " + ", ".join(sorted(pages)))


if __name__ == "__main__":
    main()
