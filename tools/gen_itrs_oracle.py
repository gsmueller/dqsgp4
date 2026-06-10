#!/usr/bin/env python
"""ERFA oracle generator for the Earth-fixed (ITRS) chain (R4b).

Emits born-digital reference values for tests/test_nutation (gate NUT1):
  - nutation angles dpsi, deps        <- erfa.nut06a   (IAU 2000A+2006, the FULL model
                                          incl. the planetary terms we bound, not model)
  - GMST / GAST (IAU 2006/2000A)      <- erfa.gmst06 / erfa.gst06a
  - polar motion matrix               <- erfa.pom00 (sp from erfa.sp00)
  - the full GCRS->ITRS matrix        <- erfa.c2t06a

Epochs are chosen on dates present in the in-repo CSSI EOP file
(sgp4_references/.../EOP-All-v1.1_2025-01-10.txt); xp, yp, dUT1, dAT are PARSED
from it (born-digital IERS data), and TT/UT1 are formed manually:
  TT  = UTC + (dAT + 32.184 s)/86400        UT1 = UTC + dUT1/86400.

Also arbitrates the GMST06 polynomial transcription: the in-repo AstroLib form
(gstime00, t^4/t^5 signs +/+) vs erfa.gmst06 — measured, recorded in the theory
note. Offline; erfa 2.0.1.5 local. Run: python tools/gen_itrs_oracle.py
"""
import math
import erfa

EOP = "sgp4_references/vallado_celestrak/datalib/EOP-All-v1.1_2025-01-10.txt"
AS2R = math.pi / (180.0 * 3600.0)

# Parse the CSSI EOP file: data lines are "yyyy mm dd mjd xp yp dut1 lod ... dat".
eop = {}
with open(EOP) as f:
    in_data = False
    for ln in f:
        s = ln.split()
        if not s:
            continue
        if s[0] in ("BEGIN", "END", "NUM_OBSERVED_POINTS", "NUM_PREDICTED_POINTS"):
            in_data = s[0] == "BEGIN"
            continue
        if not in_data or len(s) < 13:
            continue
        try:
            y, mo, d = int(s[0]), int(s[1]), int(s[2])
        except ValueError:
            continue
        eop[(y, mo, d)] = {"xp": float(s[4]), "yp": float(s[5]),
                           "dut1": float(s[6]), "dat": int(float(s[12]))}

DATES = [(2000, 1, 1), (2010, 9, 15), (2020, 1, 1), (2024, 6, 1), (1995, 3, 20)]

def jd_utc_0h(y, m, d):
    # Standard Fliegel-Van Flandern integer JD at 0h, minus 0.5 to land on 0h UTC.
    a = (14 - m) // 12
    yy = y + 4800 - a
    mm = m + 12 * a - 3
    jdn = d + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045
    return jdn - 0.5

print("# {jd_tt, jd_ut1, xp_as, yp_as, dpsi_rad, deps_rad, gmst_rad, gast_rad,")
print("#  c2t row-major 9, pom row-major 9}   (erfa nut06a/gmst06/gst06a/pom00/c2t06a)")
worst_gmst_poly = 0.0
for (y, m, d) in DATES:
    e = eop[(y, m, d)]
    jd_utc = jd_utc_0h(y, m, d)
    jd_tt = jd_utc + (e["dat"] + 32.184) / 86400.0
    jd_ut1 = jd_utc + e["dut1"] / 86400.0

    dpsi, deps = erfa.nut06a(jd_tt, 0.0)
    gmst = erfa.gmst06(jd_ut1, 0.0, jd_tt, 0.0)
    gast = erfa.gst06a(jd_ut1, 0.0, jd_tt, 0.0)
    sp = erfa.sp00(jd_tt, 0.0)
    pom = erfa.pom00(e["xp"] * AS2R, e["yp"] * AS2R, sp)
    c2t = erfa.c2t06a(jd_tt, 0.0, jd_ut1, 0.0, e["xp"] * AS2R, e["yp"] * AS2R)

    c2t9 = ", ".join(f"{c2t[i][j]:.17e}" for i in range(3) for j in range(3))
    pom9 = ", ".join(f"{pom[i][j]:.17e}" for i in range(3) for j in range(3))
    print(f"  {{{jd_tt:.10f}, {jd_ut1:.10f}, {e['xp']:.7f}, {e['yp']:.7f},")
    print(f"   {dpsi:.17e}, {deps:.17e}, {gmst:.17e}, {gast:.17e},")
    print(f"   {{{c2t9}}},")
    print(f"   {{{pom9}}}}},")

    # GMST06 polynomial arbitration: AstroLib gstime00 form (+t^4, +t^5).
    t = (jd_tt - 2451545.0) / 36525.0
    era = math.fmod(2.0 * math.pi * (0.7790572732640
                                     + 1.00273781191135448 * (jd_ut1 - 2451545.0)),
                    2.0 * math.pi)
    poly = (0.014506 + 4612.156534 * t + 1.3915817 * t**2
            - 0.00000044 * t**3 + 0.000029956 * t**4 + 0.0000000368 * t**5)
    gmst_v = math.fmod(era + poly * AS2R, 2.0 * math.pi)
    diff = abs(math.fmod(gmst_v - gmst + 3 * math.pi, 2 * math.pi) - math.pi)
    worst_gmst_poly = max(worst_gmst_poly, diff)

print(f"\n# AstroLib-form GMST06 vs erfa.gmst06: worst |diff| = "
      f"{worst_gmst_poly:.3e} rad = {worst_gmst_poly / AS2R * 1e6:.3f} uas")
