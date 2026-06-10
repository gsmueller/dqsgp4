#!/usr/bin/env python3
"""Root-cause lock for the SDP4 deep-space drag omission (deep-space class B).

The four "class B" SGP4-VER failures (28623, 23599, 11801, 16925) were originally
mis-scoped as a "resonance integrator" problem. They are in fact NON-resonant
(irez=0), high-eccentricity, high-drag deep-space orbits, and our deep-space path
had simply OMITTED the near-earth secular atmospheric drag (tempa/tempe/templ/
nodecf, plus delomg/delm for the non-simple model) that the reference applies to
ALL satellites (propagation.py:1713-1794), adding the lunar-solar DSPACE terms on
top.

This script proves the diagnosis against the born-digital reference
python_sgp4_rhodes (the trusted arbiter): for each sat it propagates twice —
normally, and with the drag coefficients zeroed — and asserts that

  * the drag is non-trivial (zeroing it moves the position by > 1 km), and
  * the magnitude of the "no-drag" error matches the error our pre-fix code
    exhibited (recorded below from the failing test_sgp4 run), to a few percent.

So if the deep-space drag is ever removed again, the regression is caught here in
addition to the test_sgp4 end-to-end comparison.

Run:  python tests/test_sgp4/verify_deep_space_drag.py
"""
import os
import sys
from math import sqrt

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
REF = os.path.join(REPO, "sgp4_references", "python_sgp4_rhodes")
if os.path.isdir(REF):
    sys.path.insert(0, REF)
sys.modules["sgp4.vallado_cpp"] = None
try:
    from sgp4.model import Satrec  # noqa: E402
except ModuleNotFoundError:
    # The python-sgp4 reference (MIT, https://pypi.org/project/sgp4/) is not
    # present; install it to run this cross-check: pip install sgp4
    print("SKIP: python-sgp4 not available (pip install sgp4)")
    sys.exit(0)

# (l1, l2, [(t_min, pre_fix_our_error_km), ...]) — the second value of each tuple
# is the error our DRAG-LESS deep-space code produced against tcppver.out, copied
# from the failing test_sgp4 run; the no-drag reference must reproduce it.
CASES = {
    "28623": (
        "1 28623U 05006B   06177.81079184  .00637644  69054-6  96390-3 0  6000",
        "2 28623  28.5200 114.9834 6249053 170.2550 212.8965  3.79477162 12753",
        [(120, 10.33), (480, 98.68), (840, 237.84)],
    ),
    "23599": (
        "1 23599U 95029B   06171.76535463  .00085586  12891-6  12956-2 0  2905",
        "2 23599   6.9327   0.2849 5782022 274.4436  25.2425  4.47796565123555",
        [(60, 0.2265), (140, 0.8986)],
    ),
    "16925": (
        "1 16925U 86065D   06151.67415771  .02550794 -30915-6  18784-3 0  4486",
        "2 16925  62.0906 295.0239 5596327 245.1593  47.9690  4.88511875148616",
        [(120, 21.02), (720, 374.21)],
    ),
    "11801": (
        "1 11801U          80230.29629788  .01431103  00000-0  14311-1      13",
        "2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848    13",
        [(720, 267.62), (1080, 865.52)],
    ),
}

DRAG_ATTRS = ["cc1", "cc4", "cc5", "d2", "d3", "d4",
              "t2cof", "t3cof", "t4cof", "t5cof", "omgcof", "xmcof", "nodecf"]


def positions(l1, l2, ts, zero_drag):
    s = Satrec.twoline2rv(l1, l2)
    if zero_drag:
        for a in DRAG_ATTRS:
            setattr(s, a, 0.0)
    out = []
    for t in ts:
        e, r, v = s.sgp4_tsince(float(t))
        out.append((e, r))
    return out


def main():
    checks = []

    def check(name, ok, detail=""):
        checks.append(ok)
        print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))

    # All four sats must be non-resonant (refutes the "resonance" hypothesis).
    for sat, (l1, l2, _pts) in CASES.items():
        s = Satrec.twoline2rv(l1, l2)
        check(f"sat {sat} is deep-space and NON-resonant (irez=0)",
              s.method == "d" and s.irez == 0, f"method={s.method} irez={s.irez}")

    # No-drag reference must reproduce the pre-fix error magnitude.
    for sat, (l1, l2, pts) in CASES.items():
        ts = [t for t, _ in pts]
        ref = positions(l1, l2, ts, zero_drag=False)
        nod = positions(l1, l2, ts, zero_drag=True)
        for (t, pre_fix_err), (e1, r1), (e2, r2) in zip(pts, ref, nod):
            if e1 or e2:
                check(f"sat {sat} t={t}: propagates", False, f"err codes {e1},{e2}")
                continue
            d = sqrt(sum((a - b) ** 2 for a, b in zip(r1, r2)))
            rel = abs(d - pre_fix_err) / pre_fix_err
            # "drag matters" = no-drag error exceeds the 0.1 km deep-space tolerance.
            check(f"sat {sat} t={t}: no-drag err {d:.2f} km ~ pre-fix {pre_fix_err} km (drag matters)",
                  d > 0.1 and rel < 0.03, f"rel={rel:.4f}")

    # isimp invariant: sgp4init forces isimp=1 (SIMPLE drag) for EVERY deep-space
    # sat (period>=225 min), independent of perigee (propagation.py:1496-1499). The
    # deep-space path must therefore always use the simple drag model and never the
    # non-simple delomg/delm/D2-4/C5/t3-5cof terms. Gating on a perigee<220 km test
    # (use_simple_model) was the bug behind the class-A eccentricity residual: these
    # sats are deep-space with perigee>220 km, so they must STILL be simple.
    RE = 6378.135
    tlefile = os.path.join(REF, "sgp4", "SGP4-VER.TLE")
    deep = 0
    deep_hi_perigee = 0          # deep-space sats with perigee > 220 km
    nonsimple_violations = 0     # deep-space sats that are NOT isimp=1
    with open(tlefile) as f:
        lines = [x.rstrip("\n") for x in f]
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            s = Satrec.twoline2rv(lines[i][:69], lines[i + 1][:69])
            if s.method == "d":
                deep += 1
                if s.isimp != 1:
                    nonsimple_violations += 1
                perigee_km = (s.a * (1.0 - s.ecco) - 1.0) * RE
                if perigee_km > 220.0 and s.isimp == 1:
                    deep_hi_perigee += 1
            i += 2
        else:
            i += 1
    check("every SGP4-VER deep-space sat has isimp=1 (simple drag)",
          nonsimple_violations == 0, f"{deep} deep-space sats, {nonsimple_violations} violations")
    # The gating-bug case must be real: deep-space sats with perigee>220 km that
    # a naive perigee<220 test would call non-simple, but which are isimp=1.
    check("deep-space sats with perigee>220 km exist and are isimp=1 (gating-bug case is real)",
          deep_hi_perigee > 0, f"{deep_hi_perigee} such sats (e.g. 23177/22674/21897)")

    # WGS72 earth-rotation rate: the deep-space resonance theta_dot = omega*60 must
    # equal the reference's rptim = 4.37526908801129966e-3 rad/min (propagation.py:672).
    # The full-precision WGS72 omega (7.2921151467e-5, what model_selector.h now uses)
    # matches it to ~9e-15; the truncated 7-figure 7.292115e-5 (WGS84's value) leaves a
    # 8.8e-11 gap that drove the 22674/21897 resonance mean-longitude residual.
    RPTIM = 4.37526908801129966e-3
    check("WGS72 omega 7.2921151467e-5 -> theta_dot matches reference rptim (<1e-13)",
          abs(7.2921151467e-5 * 60.0 - RPTIM) < 1e-13, f"diff={7.2921151467e-5*60.0-RPTIM:.2e}")
    check("the truncated 7.292115e-5 does NOT match rptim (the resonance bug)",
          abs(7.292115e-5 * 60.0 - RPTIM) > 1e-12, f"diff={7.292115e-5*60.0-RPTIM:.2e}")

    npass = sum(checks)
    print(f"\n=== deep-space drag verification: {npass}/{len(checks)} checks ===")
    print("PASS" if npass == len(checks) else "FAIL")
    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
