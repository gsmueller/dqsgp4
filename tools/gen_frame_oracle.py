#!/usr/bin/env python
"""Frame-chain oracle generator (L3 ecliptic-of-date -> GCRS).

Two independent oracle levels for test_frame_chain (design/derivations/frame_chain.md §6):

  1. ELEMENT-WISE (conformance): the IAU 2006 bias-precession matrix erfa.pmat06, the
     ecliptic<->ICRS matrix erfa.ecm06, and the mean obliquity erfa.obl06, at the four §A
     epochs. These ARE the IAU analytical theory, so matching them proves our rotations
     implement IAU 2006 (NOT a fidelity claim).
  2. END-TO-END (fidelity): the Moon's GCRS (J2000 mean-equatorial) position straight from
     the in-repo JPL DE430 table (numerical integration, columns 10-12 km) -- the TRULY
     independent ephemeris. The gate forms moon_meeus -> ecliptic Cartesian -> ecliptic_to_gcrs
     and compares the DIRECTION; the residual is the real Meeus model fidelity in GCRS.

Offline; no network. Run: `python tools/gen_frame_oracle.py`.
"""
import numpy as np
import erfa
from astropy.time import Time

# §A element-wise epochs (JD TT): 1900, J2000, ~2020.3, ~2049.8.
EPOCHS = [2415020.0, 2451545.0, 2459000.5, 2469807.5]

print("# ---- element-wise erfa oracle (pmat06, ecm06, obl06) ----")
for jd in EPOCHS:
    P = erfa.pmat06(jd, 0.0)
    M = erfa.ecm06(jd, 0.0)
    eps = erfa.obl06(jd, 0.0)
    print(f"# JD {jd:.1f}   obl06 = {eps:.18e} rad")
    print(f"  {{ {jd:.1f}, {eps:.18e},")
    print("    {" + ", ".join(f"{v:.18e}" for v in P.flatten()) + "},")
    print("    {" + ", ".join(f"{v:.18e}" for v in M.flatten()) + "} },")

# ---- end-to-end DE430 Moon GCRS (J2000 equatorial, km) ----
DE = "sgp4_references/vallado_celestrak/datalib/sunmooneph_430t12.txt"
tab = {}
with open(DE) as f:
    for ln in f:
        p = ln.split()
        if len(p) < 12:
            continue
        tab[(int(p[0]), int(p[1]), int(p[2]), int(p[3]))] = (
            float(p[9]), float(p[10]), float(p[11]))

SAMPLES = [(2000, 1, 1, 12), (2010, 9, 15, 12), (2020, 1, 1, 0),
           (2024, 6, 1, 0), (1995, 3, 20, 0)]
print("\n# ---- end-to-end DE430 Moon GCRS (JD TT, x, y, z km) ----")
for (y, mo, d, h) in SAMPLES:
    x, yk, z = tab[(y, mo, d, h)]
    jd = Time(f"{y}-{mo:02d}-{d:02d} {h:02d}:00:00", scale="tt").jd
    print(f"  {{ {jd:.6f}, {x:.4f}, {yk:.4f}, {z:.4f} }},")
