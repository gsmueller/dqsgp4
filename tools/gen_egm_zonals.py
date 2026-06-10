#!/usr/bin/env python
"""EGM2008 zonal-coefficient transcription generator (R5/O2, closes issue R04).

Reads the in-repo born-digital NGA EGM2008 normalized coefficient file
(sgp4_references/vallado_celestrak/datalib/EGM-08norm100.txt — 15-significant-
digit C_bar_n0 values) and emits, for n = 2..9:

  - the full-precision normalized C_bar_n0 strings (for
    constants/zonal_harmonics.h's generative J_n = -sqrt(2n+1)*C_bar_n0), and
  - the derived unnormalized J_n at full precision (for the
    sgp4/model_selector.h modern presets, which store J_n directly),

plus the AUDIT against the previously stored 7-digit values — the R04 record:
the old n = 5..9 transcriptions were wrong at the 5th-7th significant figure.

The gen_lunar_terms/gen_nutation_terms pattern: transcription is mechanical and
reproducible, never typed from memory. Run: python tools/gen_egm_zonals.py
"""
import math
from decimal import Decimal, getcontext

getcontext().prec = 30

EGM = "sgp4_references/vallado_celestrak/datalib/EGM-08norm100.txt"

# The previously stored values (the R04 state), for the audit record.
OLD_CBAR = {3: "0.9571612e-6", 5: "0.0686729e-6", 7: "0.0905123e-6", 9: "0.0280234e-6"}
OLD_JN = {3: "-0.00000253215306", 4: "-0.00000161989770", 5: "-0.000000227735",
          6: "0.000000540774", 7: "-0.000000350447", 8: "-0.000000203948",
          9: "-0.000000122154"}

cbar = {}
with open(EGM) as f:
    for ln in f:
        p = ln.split()
        if len(p) < 4:
            continue
        n, m = int(p[0]), int(p[1])
        if m == 0 and 2 <= n <= 9:
            cbar[n] = p[2]

print("# n   C_bar_n0 (file, verbatim)        J_n = -sqrt(2n+1)*C_bar_n0 (15 digits)")
for n in range(2, 10):
    c = Decimal(cbar[n])
    jn = -Decimal(math.isqrt((2 * n + 1) * 10**40)).scaleb(-20) * c  # sqrt via isqrt, 20 digits
    print(f"  {n}   {cbar[n]:<28}   {jn:.15e}")

print("\n# R04 audit: old stored values vs the file")
for n in sorted(OLD_CBAR):
    old = Decimal(OLD_CBAR[n])
    new = Decimal(cbar[n])
    rel = abs(old - new) / abs(new)
    print(f"  Cbar{n}0: stored {OLD_CBAR[n]:<14} file {cbar[n]:<24} rel diff {float(rel):.2e}")
for n in sorted(OLD_JN):
    c = Decimal(cbar[n])
    jn = -Decimal(math.isqrt((2 * n + 1) * 10**40)).scaleb(-20) * c
    old = Decimal(OLD_JN[n])
    rel = abs(old - jn) / abs(jn)
    print(f"  J{n}: stored {OLD_JN[n]:<18} derived {jn:.9e}  rel diff {float(rel):.2e}")
