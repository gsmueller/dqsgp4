#!/usr/bin/env python
"""JPL DE430 third-body oracle generator (L4, INDEPENDENT of Meeus).

Emits the Sun's and Moon's geocentric positions (J2000 mean-equatorial ≈ GCRS, km)
from the in-repo JPL DE430 ephemeris table — the numerical-integration reference the
Cartesian third-body force (src/forces/third_body.h) is gated against in
tests/test_third_body. The force feeds the MODEL ephemeris (sun_meeus / moon_meeus)
through the same Battin acceleration; comparing to the acceleration from these DE430
positions isolates the ephemeris-truncation residual (the body-position accuracy must
majorize it — no perceived fidelity). DE430 is a numerical integration, so it is truly
independent of the analytical Meeus/ELP theory the model realises (the independence
lesson — astropy's builtin Sun/Moon are erfa.epv00/moon98, the SAME theory).

The table (sgp4_references/vallado_celestrak/datalib/sunmooneph_430t12.txt, 12-hour
steps 1956–) holds geocentric J2000-equatorial positions in km: Sun xyz in columns
5–7 (0-indexed 4–6), Moon xyz in columns 10–12 (0-indexed 9–11). Offline; no network.
Run: `python tools/gen_third_body_oracle.py`.
"""
from astropy.time import Time

DE = "sgp4_references/vallado_celestrak/datalib/sunmooneph_430t12.txt"
tab = {}
with open(DE) as f:
    for ln in f:
        p = ln.split()
        if len(p) < 12:
            continue
        tab[(int(p[0]), int(p[1]), int(p[2]), int(p[3]))] = (
            float(p[4]), float(p[5]), float(p[6]),     # Sun xyz [km]
            float(p[9]), float(p[10]), float(p[11]))   # Moon xyz [km]

# Sample epochs present in the 12-hour grid, spanning decades.
SAMPLES = [(2000, 1, 1, 12), (2010, 9, 15, 12), (2020, 1, 1, 0), (2024, 6, 1, 0), (1995, 3, 20, 0)]

print("# {JD(TT), sunX,sunY,sunZ [km], moonX,moonY,moonZ [km]}  (J2000 mean-equatorial ~ GCRS)")
for (y, mo, d, h) in SAMPLES:
    sx, sy, sz, mx, my, mz = tab[(y, mo, d, h)]
    jd = Time(f"{y}-{mo:02d}-{d:02d} {h:02d}:00:00", scale="tt").jd
    print(f"  {{{jd:.1f}, {sx:.4f},{sy:.4f},{sz:.4f}, {mx:.4f},{my:.4f},{mz:.4f}}},")
