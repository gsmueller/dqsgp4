#!/usr/bin/env python
"""Astropy solar-ephemeris oracle generator (L3, no perceived fidelity).

Emits the Sun's GEOMETRIC geocentric ecliptic-of-date longitude and distance at
sample epochs, as the INDEPENDENT reference that src/ephemeris/sun_meeus.h is
gated against (test_ephemeris / EPH). astropy/ERFA is independent of any in-repo
code, so a match is proof — the oracle is established FIRST, before any fidelity
is claimed (design/derivations/ephemeris.md §6).

astropy `get_body` returns the LIGHT-TIME-CORRECTED apparent position; the model
(and the third-body perturbation) want the instantaneous GEOMETRIC direction, so
we undo the annual aberration with the IAU constant kappa = 20.4955" (kappa/R in
longitude). The residual after this correction is the true Meeus §25 truncation
grade (<= ~0.007 deg). Run: `python tools/gen_solar_oracle.py`. Offline (builtin
ERFA/VSOP ephemeris). Paste the printed rows into test_ephemeris's oracle table.
"""
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body, GeocentricMeanEcliptic, solar_system_ephemeris

solar_system_ephemeris.set('builtin')
KAPPA_ARCSEC = 20.4955  # IAU constant of aberration

# Sample epochs (JD, TT): J2000, and four 21st-century dates spanning the orbit.
JDS = [2451545.0, 2455000.5, 2459000.5, 2460310.5, 2462000.5]

print(f"# astropy {__import__('astropy').__version__}, builtin ephemeris; geometric = apparent + kappa/R")
print(f"# {'JD (TT)':>12}  {'lon_geom [deg]':>16}  {'R [AU]':>12}")
for jd in JDS:
    t = Time(jd, format='jd', scale='tt')
    ecl = get_body('sun', t).transform_to(GeocentricMeanEcliptic(equinox=t))
    R = ecl.distance.to('AU').value
    lon_geom = (ecl.lon.deg + (KAPPA_ARCSEC / R) / 3600.0) % 360.0
    print(f"  {jd:12.1f}  {lon_geom:16.9f}  {R:12.9f}")
