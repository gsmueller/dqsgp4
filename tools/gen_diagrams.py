#!/usr/bin/env python
"""Generated SVG diagrams for the documentation set (Q3, 2026-06-10).

Imported by gen_docs.py; each diagram is embedded inline in the pages, so the
DOC1 regenerate-and-diff gate covers diagram freshness automatically.

Honesty rules:
  - The layer/dependency map's EDGES are scanned from the real #include graph
    at generation time (they cannot rot); only the column layout is curated.
  - Every numeric grade on an annotation is a MEASURED gate result (the gate
    is named beside it); pipelines and groupings are structural, not graded.

Standalone preview: python tools/gen_diagrams.py  (writes build/diagrams_preview.html)
"""
import html
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")

# The docs palette (matches gen_docs CSS).
INK = "#1a232e"
DEEP = "#102a43"
MID = "#243b53"
ACCENT = "#486581"
LINE = "#9fb3c8"
PANE = "#f0f4f8"
EDGE = "#d9e2ec"
WARN = "#8d6708"
GOOD = "#1f7a4d"


def _esc(s):
    return html.escape(str(s), quote=True)


class Svg:
    """Tiny SVG builder: boxes with multi-line labels + marker arrows."""

    def __init__(self, name, w, h):
        self.name = name
        self.w = w
        self.h = h
        self.parts = []

    def box(self, x, y, w, h, lines, fill=PANE, stroke=ACCENT, size=13,
            bold_first=True, text_fill=INK):
        self.parts.append(
            f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="7"'
            f' fill="{fill}" stroke="{stroke}" stroke-width="1.4"/>')
        if isinstance(lines, str):
            lines = [lines]
        n = len(lines)
        for i, ln in enumerate(lines):
            ty = y + h / 2 + (i - (n - 1) / 2) * (size + 3) + size * 0.35
            weight = ' font-weight="600"' if (i == 0 and bold_first) else ""
            fsize = size if i == 0 else size - 2
            self.parts.append(
                f'<text x="{x + w / 2}" y="{ty:.1f}" text-anchor="middle"'
                f' font-size="{fsize}"{weight} fill="{text_fill}">{_esc(ln)}</text>')

    def label(self, x, y, text, size=11, fill=ACCENT, anchor="middle", bold=False):
        weight = ' font-weight="600"' if bold else ""
        self.parts.append(
            f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}"'
            f'{weight} fill="{fill}">{_esc(text)}</text>')

    def arrow(self, x1, y1, x2, y2, color=ACCENT, dashed=False, width=1.6,
              both=False):
        dash = ' stroke-dasharray="6 4"' if dashed else ""
        m0 = f' marker-start="url(#{self.name}-arr)"' if both else ""
        self.parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}"'
            f' stroke="{color}" stroke-width="{width}"{dash}'
            f' marker-end="url(#{self.name}-arr)"{m0}/>')

    def render(self):
        return (
            f'<svg viewBox="0 0 {self.w} {self.h}" role="img"'
            f' xmlns="http://www.w3.org/2000/svg"'
            f' style="max-width:{self.w}px;width:100%;height:auto;display:block;'
            f'margin:14px auto;font-family:Segoe UI,Helvetica,Arial,sans-serif">'
            f'<defs><marker id="{self.name}-arr" viewBox="0 0 10 10" refX="9"'
            f' refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M 0 0 L 10 5 L 0 10 z" fill="{ACCENT}"/></marker></defs>'
            + "".join(self.parts) + "</svg>")


# ---------------------------------------------------------------- layer map

def include_graph():
    """module -> set(directly included modules), scanned from src/ itself."""
    mods = [d for d in os.listdir(SRC) if os.path.isdir(os.path.join(SRC, d))]
    edges = {}
    for m in mods:
        deps = set()
        for fn in os.listdir(os.path.join(SRC, m)):
            if not fn.endswith((".h", ".cpp")):
                continue
            try:
                with open(os.path.join(SRC, m, fn), encoding="utf-8",
                          errors="replace") as f:
                    text = f.read()
            except OSError:
                continue
            for dep in re.findall(r'#include\s+"\.\./([a-z_0-9]+)/', text):
                if dep != m and dep in mods:
                    deps.add(dep)
        edges[m] = deps
    return edges


# Curated COLUMNS only (foundation -> surface); the edges are scanned live.
LAYER_COLUMNS = [
    ["math"],
    ["geodesy", "tle", "perturbation"],
    ["constants", "astronomy"],
    ["ephemeris", "orbit", "atmosphere"],
    ["forces", "integrators", "sgp4"],
    ["dynamics"],
]


def layer_map():
    edges = include_graph()
    col_of = {m: ci for ci, col in enumerate(LAYER_COLUMNS) for m in col}
    col_w, box_w, box_h = 168, 138, 46
    rows = max(len(c) for c in LAYER_COLUMNS)
    h = 96 + rows * 86
    s = Svg("layers", 40 + col_w * len(LAYER_COLUMNS), h)
    pos = {}
    for ci, col in enumerate(LAYER_COLUMNS):
        x = 24 + ci * col_w
        for ri, m in enumerate(col):
            y = 54 + ri * 86 + (rows - len(col)) * 43 // 1
            pos[m] = (x, y)
    s.label(s.w / 2, 22, "Module layer map — edges scanned from the real "
            "#include graph at generation time", size=13, fill=DEEP, bold=True)
    s.label(s.w / 2, 40, "adjacent-layer dependencies drawn (arrow points at "
            "the dependency; dashed amber = mutual header-level pair); the "
            "complete edge list is tabulated below", size=11)
    drawn = set()
    for m, deps in sorted(edges.items()):
        for d in sorted(deps):
            if d == "math" or m not in pos or d not in pos:
                continue
            # Only adjacent columns (or the same column): these arrows live in
            # the inter-column gutter by construction, so they cannot pass
            # through any box. The table below carries the full edge list.
            if abs(col_of[m] - col_of[d]) > 1:
                continue
            mutual = m in edges.get(d, set())
            key = tuple(sorted((m, d)))
            if mutual and key in drawn:
                continue
            drawn.add(key)
            x1, y1 = pos[m]
            x2, y2 = pos[d]
            sx, sy = x1 + (0 if x2 < x1 else box_w), y1 + box_h / 2
            ex, ey = x2 + (box_w if x2 < x1 else 0), y2 + box_h / 2
            if x1 == x2:
                sx = ex = x1 + box_w / 2
                sy = y1 + (box_h if y2 > y1 else 0)
                ey = y2 + (0 if y2 > y1 else box_h)
            s.arrow(sx, sy, ex, ey, color=LINE if not mutual else WARN,
                    dashed=mutual, width=1.3, both=mutual)
    for m, (x, y) in pos.items():
        nh = len([f for f in os.listdir(os.path.join(SRC, m))
                  if f.endswith(".h")])
        s.box(x, y, box_w, box_h, [m, f"{nh} headers"],
              fill=PANE if m != "math" else "#e3ecf3",
              stroke=DEEP if m in ("math", "dynamics", "sgp4") else ACCENT)
    s.label(s.w / 2, h - 14, "foundation → surface (left to right); "
            "the dqsgp4.h umbrella re-exports the public face of every column",
            size=11)

    # The COMPLETE scanned dependency list (mechanical), plus the mutual pairs.
    mutual_pairs = sorted({tuple(sorted((m, d)))
                           for m, deps in edges.items() for d in deps
                           if m in edges.get(d, set())})
    rows_html = "".join(
        f"<tr><td><code>{m}</code></td>"
        f"<td>{', '.join(f'<code>{d}</code>' for d in sorted(deps)) or '—'}</td></tr>"
        for m, deps in sorted(edges.items()))
    table = ("<table><tr><th>Module</th><th>directly includes (scanned)</th></tr>"
             + rows_html + "</table>"
             "<p class='muted'>Header-level mutual pairs (templates, "
             "header-only): "
             + "; ".join(f"<code>{a}</code> ↔ <code>{b}</code>"
                         for a, b in mutual_pairs)
             + ".</p>")
    return s.render() + table


# ------------------------------------------------------------- frame chain

def frame_chain():
    s = Svg("frames", 1040, 300)
    s.label(520, 24, "The frame chain — from the propagators' TEME to the "
            "Earth-fixed ITRS", size=14, fill=DEEP, bold=True)
    bw, bh, y = 150, 54, 60
    xs = [30, 250, 470, 690, 880]
    s.box(xs[0], y, bw, bh, ["TEME", "SGP4 / DQ output"])
    s.box(xs[1], y, bw, bh, ["GCRS", "quasi-inertial"])
    s.box(xs[2], y, bw, bh, ["true equator,", "true equinox of date"])
    s.box(xs[3], y, bw, bh, ["ITRS", "Earth-fixed"])
    s.box(xs[4], y, 130, bh, ["subpoint", "lat / lon / alt"])
    s.arrow(xs[0] + bw, y + bh / 2, xs[1], y + bh / 2)
    s.arrow(xs[1] + bw, y + bh / 2, xs[2], y + bh / 2)
    s.arrow(xs[2] + bw, y + bh / 2, xs[3], y + bh / 2)
    s.arrow(xs[3] + bw, y + bh / 2, xs[4], y + bh / 2)
    s.label((xs[0] + bw + xs[1]) / 2, y - 12, "Pᵀ (inverse IAU2006 precession)")
    s.label((xs[0] + bw + xs[1]) / 2, y + bh + 18,
            "omitted nutation+EE ≤ ~24″ — deposited as an", size=10)
    s.label((xs[0] + bw + xs[1]) / 2, y + bh + 31,
            "accuracy bound (FRAME2 erfa-gated rotations)", size=10)
    s.label((xs[1] + bw + xs[2]) / 2, y - 12, "N·P (IAU 2000A nutation)")
    s.label((xs[1] + bw + xs[2]) / 2, y + bh + 18,
            "678-term in-repo SOFA table; vs erfa nut06a", size=10)
    s.label((xs[1] + bw + xs[2]) / 2, y + bh + 31,
            "0.48/0.38 mas (NUT1, planetary floor)", size=10)
    s.label((xs[2] + bw + xs[3]) / 2, y - 12, "R3(GAST) · then W(xp, yp)")
    s.label((xs[2] + bw + xs[3]) / 2, y + bh + 18,
            "GAST 0.44 mas; polar motion BIT-exact vs", size=10)
    s.label((xs[2] + bw + xs[3]) / 2, y + bh + 31,
            "erfa pom00 (NUT1)", size=10)
    s.label((xs[3] + bw + xs[4]) / 2, y - 12, "geocentric")
    s.label((xs[3] + bw + xs[4]) / 2, y + bh + 18, "EX1 ground track", size=10)
    s.box(30, 200, 980, 58, [
        "Full chain vs erfa c2t06a: 1.85e-9 per element ≈ 0.4 mas (gate NUT1, "
        "planetary-floor-dominated).",
        "The frozen SGP4 path keeps its historical Aoki-82 TEME→ECEF = "
        "Rz(GMST) convention (untouched, OR1)."],
        fill="#eef2f6", stroke=EDGE, size=12, bold_first=False)
    return s.render()


# -------------------------------------------------------- force composition

def force_composition():
    s = Svg("forces", 1040, 430)
    s.label(520, 24, "Force composition — opt-in members summing into the DQ "
            "propagator", size=14, fill=DEEP, bold=True)
    fy, fh, fw = 52, 58, 232
    s.box(30, fy, fw, fh, ["geopotential", "monopole+zonal+tesseral",
                           "one Cunningham V/W pass"], stroke=DEEP)
    s.label(30 + fw / 2, fy + fh + 14,
            "GEOPOT: vs legacy sum 4.9e-16, closed-form J₂ 2.5e-16", size=10)
    s.box(30, fy + 100, fw, fh, ["drag", "½ρB|v_rel|v_rel, pluggable "
                                 "DensityModel"])
    s.label(30 + fw / 2, fy + 100 + fh + 14,
            "density values verified against the published table", size=10)
    s.box(30, fy + 200, fw, fh, ["third body (Sun, Moon)", "Battin f(q), "
                                 "cancellation-free"])
    s.label(30 + fw / 2, fy + 200 + fh + 14,
            "TB1: ephemeris→accel vs JPL DE430", size=10)
    s.box(30, fy + 300, fw, fh, ["solar radiation pressure", "cannonball + "
                                 "cylindrical shadow"])
    s.label(30 + fw / 2, fy + 300 + fh + 14,
            "SRP1: P₁ᴬᵁ GENERATED from L☉/(4π·AU²·c)", size=10)

    s.box(360, 150, 150, 64, ["Σ Wrench", "force list"], fill="#e3ecf3",
          stroke=DEEP)
    for yy in (fy + fh / 2, fy + 100 + fh / 2, fy + 200 + fh / 2, fy + 300 + fh / 2):
        s.arrow(30 + fw, yy, 360, 182)
    s.box(570, 150, 190, 64, ["integrator", "RK4 fixed / RKF7(8) adaptive"])
    s.arrow(510, 182, 570, 182)
    s.label(665, 230, "RK1 order conditions exact;", size=10)
    s.label(665, 243, "AD1 adaptive bit-identical to the standalone loop", size=10)
    s.box(820, 150, 190, 64, ["State<T>", "SE(3) dual-quaternion"])
    s.arrow(760, 182, 820, 182)

    s.box(360, 300, 650, 100, [
        "Default model = geopotential only (bit-frozen, FM1). Presets compose: "
        "lunisolar / drag_B / SRP via DqForceOptions;",
        "arbitrary extra_forces append after the core. FM1 measures the presets "
        "on satellite 00005 over 60 min:",
        "lunisolar 12.6 m  >  SRP 1.39 m  >  drag 0.33 m — and asserts the "
        "default stays bit-identical."],
        fill="#eef2f6", stroke=EDGE, size=12, bold_first=False)
    return s.render()


# --------------------------------------------------------- three-error flow

def three_error_flow():
    s = Svg("errors", 1040, 360)
    s.label(520, 24, "The three-error budget through one operation",
            size=14, fill=DEEP, bold=True)
    s.box(30, 60, 250, 110, ["TrackedValue<T> a", "value",
                             "measurement σ", "precision δp",
                             "accuracy δa"], size=13)
    s.box(30, 200, 250, 110, ["TrackedValue<T> b", "value",
                              "measurement σ", "precision δp",
                              "accuracy δa"], size=13)
    s.box(380, 130, 220, 110, ["operation", "(+, ×, sin, exp, …)"],
          fill="#e3ecf3", stroke=DEEP)
    s.arrow(280, 115, 380, 170)
    s.arrow(280, 255, 380, 200)
    s.box(740, 130, 270, 110, ["result", "value: the arithmetic",
                               "each channel: first-order propagation",
                               "+ this op's own rounding → precision"],
          size=13)
    s.arrow(600, 185, 740, 185)
    s.box(30, 312, 980, 44, [
        "measurement: physical input σ only · precision: representation/rounding "
        "(tightens with wider T)",
        "accuracy: model truncation via add_bound (never improved by arithmetic) "
        "· total() = the conservative sum"],
        fill="#eef2f6", stroke=EDGE, size=11, bold_first=False)
    s.label(520, 295, "Transcendentals carry sup-derivative bounds; division "
            "with an error interval containing zero yields the honest "
            "max-sentinel; gate W12 asserts the propagation.", size=11)
    return s.render()


# ------------------------------------------------------- atmosphere ladder

def atmosphere_ladder():
    s = Svg("atmo", 1040, 330)
    s.label(520, 24, "The atmosphere model ladder — what each rung adds",
            size=14, fill=DEEP, bold=True)
    bw, bh = 300, 92
    s.box(30, 60, bw, bh, ["SGP4 Lane power-law", "frozen inside the analytical "
                           "path", "(s, q-s)⁴ profile; OR1 bit-frozen"],
          size=13)
    s.box(370, 60, bw, bh, ["Vallado 8-4 static table", "27 exponential bands, "
                            "born-digital", "ATM1: nodes bit-exact; chaining "
                            "9.6e-5"], size=13)
    s.box(710, 60, bw, bh, ["space-weather model", "(e.g. NRLMSISE-00) plugs the",
                            "same DensityModel interface"], size=13, stroke=DEEP)
    s.arrow(30 + bw, 106, 370, 106)
    s.arrow(370 + bw, 106, 710, 106)
    s.label(350, 86, "absolute altitude", size=10)
    s.label(350 + 340, 86, "solar activity", size=10)
    s.box(30, 190, 980, 110, [
        "What a static table cannot represent (measured against an "
        "NRLMSISE-00 reference, 2026-06-09):",
        "solar-cycle density range at 400 km = 11.7× (71× at 700 km); "
        "diurnal bulge 2.5×; T∞ 884/1251/1604 K low/moderate/high.",
        "A space-weather model takes the F10.7a / F10.7 / Ap indices and "
        "carries its published density error as its accuracy band.",
        "The DensityModel interface means every drag consumer picks the rung "
        "without code changes (write the model, pass it to make_drag)."],
        fill="#eef2f6", stroke=EDGE, size=12, bold_first=False)
    return s.render()


# ------------------------------------------------------ sgp4 vs dqsgp4 split

def sgp4_vs_dqsgp4():
    s = Svg("split", 1040, 430)
    s.label(520, 24, "From one TLE to either propagator",
            size=14, fill=DEEP, bold=True)
    s.box(420, 50, 200, 50, ["TLE / OMM", "tle::parse"], fill="#e3ecf3",
          stroke=DEEP)
    # left: SGP4
    s.box(80, 150, 360, 96, ["sgp4::Propagator<T>", "authentic analytical "
                             "SGP4/SDP4 (WGS72)", "TEME kilometres"],
          stroke=DEEP)
    s.label(260, 262, "verified: 33/33 satellites, 623/623 points of the",
            size=11, fill=GOOD)
    s.label(260, 276, "published SGP4 test set; results pinned bit-for-bit", size=11,
            fill=GOOD)
    # right: DQSGP4
    s.box(600, 150, 360, 96, ["dynamics::DqSgp4Propagator<T>",
                              "numerical SE(3) integration, metres",
                              "authentic | boosted modes + presets"],
          stroke=DEEP)
    s.label(780, 262, "initial state recovered through the WGS72 SGP4 model",
            size=11)
    s.label(780, 276, "(a TLE's elements are defined by it); forces opt-in", size=11)
    s.arrow(470, 100, 260, 150)
    s.arrow(570, 100, 780, 150)
    s.arrow(440, 198, 600, 198, dashed=True)
    s.label(520, 190, "seed", size=10)
    # shared verb
    s.box(330, 320, 380, 64, ["the same verb",
                              "propagate(minutes since epoch)",
                              "Propagatable concept"], fill="#e3ecf3",
          stroke=DEEP)
    s.arrow(260, 246, 430, 320)
    s.arrow(780, 246, 610, 320)
    s.label(520, 404, "The test suite drives both propagators through the shared "
            "interface and verifies the results match direct calls exactly.",
            size=11)
    return s.render()


def diagrams():
    """name -> inline-embeddable SVG markup."""
    return {
        "layer_map": layer_map(),
        "frame_chain": frame_chain(),
        "force_composition": force_composition(),
        "three_error_flow": three_error_flow(),
        "atmosphere_ladder": atmosphere_ladder(),
        "sgp4_vs_dqsgp4": sgp4_vs_dqsgp4(),
    }


if __name__ == "__main__":
    body = "".join(f"<h2>{n}</h2>{svg}" for n, svg in diagrams().items())
    out = os.path.join(ROOT, "build", "diagrams_preview.html")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        f.write('<!DOCTYPE html><html><head><meta charset="utf-8"></head>'
                f"<body>{body}</body></html>")
    print(f"wrote {out}")
