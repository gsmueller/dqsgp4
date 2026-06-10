#!/usr/bin/env python3
"""Reference-behavior lock for the SDP4 error-code edge cases (deep-space class C).

The C++ test harness (tests/test_sgp4/main.cpp) must agree with the *born-digital*
reference python_sgp4_rhodes on WHERE the propagator errors. This script runs that
trusted reference over the full SGP4-VER set and asserts:

  1. the canonical error-code sequence  [1, 1, 6, 6, 4, 3, 6]  (the same list the
     reference's own tests.py:600/609/635 asserts), and
  2. sat 33334 errors at t=0 with code 3 (perturbed eccentricity out of range),
     perturbed e == -122.217193 (the value our deep_space.h guard must reject).

Only 33334 produces a t=0 error (a stale-buffer "(Use previous data line)" in
tcppver.out); the other six errors are t>0 and merely TRUNCATE their satellite's
reference output, so the position-comparison harness never reaches them. This is
the reference truth that the deep_space.h NaN guard + tcppver_parser.h stale-buffer
detection are matched against (Standard 10 for an error condition rather than a
coefficient: the code is the arbiter).

Run:  python tests/test_sgp4/verify_deep_space_error_codes.py
"""
import os
import sys
from math import isnan

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
REF = os.path.join(REPO, "sgp4_references", "python_sgp4_rhodes")
TLE = os.path.join(REF, "sgp4", "SGP4-VER.TLE")
if not os.path.exists(TLE):
    TLE = os.path.join(REPO, "sgp4_references", "aholinch_sgp4", "data", "SGP4-VER.TLE")
if os.path.isdir(REF):
    sys.path.insert(0, REF)

# Force the pure-python implementation (the born-digital reference), not the C ext.
sys.modules["sgp4.vallado_cpp"] = None
try:
    from sgp4.model import Satrec  # noqa: E402
except ModuleNotFoundError:
    # The python-sgp4 reference (MIT, https://pypi.org/project/sgp4/) is not
    # present; install it to run this cross-check: pip install sgp4
    print("SKIP: python-sgp4 not available (pip install sgp4)")
    sys.exit(0)

EXPECTED_ERRORS = [1, 1, 6, 6, 4, 3, 6]   # reference tests.py expected_errors


def run():
    error_list = []
    sat33334 = None
    with open(TLE) as f:
        lines = iter(f.read().splitlines())
    for l1 in lines:
        if not l1.startswith("1"):
            continue
        l2 = next(lines)
        sat = Satrec.twoline2rv(l1, l2)
        tstart, tend, tstep = (float(x) for x in l2[69:].split())
        # t=0 first (matches generate_satellite_output ordering)
        e, r, v = sat.sgp4_tsince(0.0)
        if isnan(r[0]):
            error_list.append(e)
            if str(sat.satnum) == "33334":
                sat33334 = (e, sat.error_message)
            continue
        t = tstart if tstart != 0.0 else tstart + tstep
        while t <= tend + 1e-6:
            e, r, v = sat.sgp4_tsince(t)
            if e != 0:
                error_list.append(e)
                break
            t += tstep
    return error_list, sat33334


def main():
    checks = []

    def check(name, ok, detail=""):
        checks.append(ok)
        print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))

    error_list, sat33334 = run()

    check("canonical error-code sequence == [1,1,6,6,4,3,6]",
          error_list == EXPECTED_ERRORS, f"got {error_list}")
    check("exactly 7 error events", len(error_list) == 7, f"got {len(error_list)}")
    check("exactly one error-3 (perturbed ecc) event",
          error_list.count(3) == 1, f"count={error_list.count(3)}")
    check("sat 33334 errors at t=0 with code 3",
          sat33334 is not None and sat33334[0] == 3, f"{sat33334}")
    if sat33334 is not None:
        msg = sat33334[1] or ""
        check("sat 33334 perturbed e == -122.217193",
              "-122.217193" in msg, msg)

    npass = sum(checks)
    print(f"\n=== deep-space error-code verification: {npass}/{len(checks)} checks ===")
    print("PASS" if npass == len(checks) else "FAIL")
    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
