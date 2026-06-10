#!/usr/bin/env python
"""Static-atmosphere table transcription generator (R1, design/derivations/
atmosphere_exponential_table.md section 3).

Reads the in-repo born-digital Vallado Table 8-4 data file (ATMOSEXP.DAT: base
altitude [km], nominal density [kg/m^3], scale height [km], 27 rows) and EMITS the
C++ rows for src/atmosphere/exponential_table.h — the gen_lunar_terms.py pattern:
transcription is mechanical and reproducible, never typed from memory.

Also reports the NEIGHBOR-CHAINING residuals (row i's exponential evaluated at row
i+1's base vs row i+1's rho0) — the table was constructed to chain, so the max
residual calibrates the gate tolerance and the model's default declared band.

Offline; no network. Run: `python tools/gen_atmosphere_table.py`.
"""
import math

DAT = "sgp4_references/vallado_celestrak/software/misc/pascal/ATMOSEXP.DAT"

rows = []
with open(DAT) as f:
    for ln in f:
        p = ln.split()
        if len(p) != 3:
            continue
        h0_km, rho0, h_km = float(p[0]), float(p[1]), float(p[2])
        rows.append((h0_km, rho0, h_km, p[1]))

print(f"# {len(rows)} rows from {DAT}")
print("# C++ rows: { h0 [m, exact int], rho0 [kg/m^3, model_coefficient], H [km, model_coefficient] }")
for (h0_km, rho0, h_km, rho0_raw) in rows:
    h0_m = round(h0_km * 1000)
    assert abs(h0_m - h0_km * 1000) < 1e-9, "base altitude not exact at the metre"
    # normalize the Pascal exponent form 3.899E-0002 -> 3.899e-2
    mant, _, expo = rho0_raw.upper().partition("E")
    rho0_str = mant + (f"e{int(expo)}" if expo and int(expo) != 0 else "")
    print(f'    {{{h0_m}, "{rho0_str}", "{h_km:.3f}"}},')

print("\n# neighbor-chaining report: rho_i(h0_{i+1}) vs rho0_{i+1}")
worst = 0.0
for i in range(len(rows) - 1):
    h0, rho0, H, _ = rows[i]
    h1, rho1, _, _ = rows[i + 1]
    chained = rho0 * math.exp(-(h1 - h0) / H)
    rel = abs(chained - rho1) / rho1
    worst = max(worst, rel)
    print(f"  {h0:7.2f} -> {h1:7.2f} km : chained {chained:.4e}  table {rho1:.4e}  rel {rel:.2e}")
print(f"# max chaining residual: {worst:.3e}")
