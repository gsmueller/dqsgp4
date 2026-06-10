#!/usr/bin/env python
"""JPL DE430 lunar oracle generator (L3, INDEPENDENT of Meeus/ELP).

Emits the Moon's geocentric ecliptic-of-date longitude/latitude/distance from the
in-repo JPL DE430 ephemeris table (numerical integration) — the reference that
src/ephemeris/moon_meeus.h is gated against (test_ephemeris). This is a TRULY
independent oracle: astropy's `builtin` Moon is itself the Meeus/ELP truncation
(erfa.moon98), so matching it only checks the code reproduces the theory; DE430
is a numerical integration, so the residual is the real model fidelity.

The DE430 table (sgp4_references/vallado_celestrak/datalib/sunmooneph_430t12.txt,
12-hour steps 1956–) gives the Moon's geocentric position in the J2000 mean
equatorial frame (km, columns 10–12). We form RA/Dec and rotate to ecliptic OF
DATE with ERFA's IAU-2006 eqec06 (the same standard our frame chain targets).
Offline; no network. Run: `python tools/gen_lunar_oracle.py`.
"""
import numpy as np
import erfa
from astropy.time import Time

DE = "sgp4_references/vallado_celestrak/datalib/sunmooneph_430t12.txt"
tab = {}
with open(DE) as f:
    for ln in f:
        p = ln.split()
        if len(p) < 12:
            continue
        tab[(int(p[0]), int(p[1]), int(p[2]), int(p[3]))] = (float(p[9]), float(p[10]), float(p[11]))

# Sample epochs present in the table (12-hour grid), spanning decades.
SAMPLES = [(2000, 1, 1, 12), (2010, 9, 15, 12), (2020, 1, 1, 0), (2024, 6, 1, 0), (1995, 3, 20, 0)]

print(f"# JPL DE430 (in-repo) -> ecliptic-of-date via erfa.eqec06.  {'JD (TT)':>14} {'lon[deg]':>13} {'lat[deg]':>12} {'dist[km]':>14}")
for (y, mo, d, h) in SAMPLES:
    x, yk, z = tab[(y, mo, d, h)]
    jd = Time(f"{y}-{mo:02d}-{d:02d} {h:02d}:00:00", scale="tt").jd
    ra = np.arctan2(yk, x)
    dec = np.arctan2(z, np.hypot(x, yk))
    lon, lat = erfa.eqec06(jd, 0.0, ra, dec)
    dist = np.sqrt(x * x + yk * yk + z * z)
    print(f"  {{{jd:.1f}, {np.degrees(lon) % 360:.7f}, {np.degrees(lat):+.7f}, {dist:.3f}}},")
