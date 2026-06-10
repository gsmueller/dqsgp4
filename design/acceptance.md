# Acceptance gates — the mechanical definition of "complete"

This document is the human-readable companion to **`tools/run_acceptance.ps1`**, which is the executable
source of truth. The project is **COMPLETE if and only if `run_acceptance.ps1` prints `COMPLETE`** (gate
`G*` in `recursive-forging-rain.md`). At any moment its verdict line — `COMPLETE` or `k gate(s)
remaining` — is the exact, objective project status. There is no other criterion; "complete" is never a
judgment call.

```
pwsh tools/run_acceptance.ps1              # full run (includes the slow Octave/Python verifiers)
pwsh tools/run_acceptance.ps1 -SkipVerifiers   # fast inner-loop run (regression + static gates only)
```

Exit code 0 ⟺ `COMPLETE`, else 1 (so the script also serves as a CI gate).

## I0 — the standing regression invariant (must always hold; preemptive)

| Gate | Check |
|---|---|
| `I0.build` | `sgp4.sln` builds Release\|x64 (full-path MSBuild), `MSBUILD_EXIT=0` |
| `I0.sgp4` | `test_sgp4.exe` prints `Satellites: 33 pass, 0 fail` and `Data points: 623 pass, 0 fail` |
| `I0.<module>` | each of test_math/geodesy/wgs84/astronomy/perturbation/tle/propagator/dual_number/quaternion/dual_quaternion exits 0 |
| `I0.verify` | every `design/**/verify_*.m` runs clean under Octave |
| `I0.pyverify` | every `tests/test_sgp4/verify_*.py` exits 0 |

If any `I0.*` is red, **restoring it is the current task** — it preempts the Work Queue until green.

## W — one gate per Work-Queue item (see `recursive-forging-rain.md` for the items)

| Gate | DoD check (what `run_acceptance` looks for) |
|---|---|
| `W1`–`W5` | `verify_near_{short_period,longperiod,kepler,secular,osculating}.m` exist and run clean (P1 fidelity) |
| `W6`,`W7` | `src/forces/gravity_central.h` and `gravity_zonal.h` write `errors.accuracy` |
| `W8` | `test_force_models.exe` exits 0 |
| `W9`,`W10` | `tools/check_index.ps1` and `tools/check_traceability.ps1` exit 0 |
| `W11`,`W12` | `test_code_consistency.exe`, `test_error_framework.exe` exit 0 |
| `W13`–`W15` | `verify_dim_{brouwer,kaula,recovery}.m` run clean |
| `W16` | `test_model_value.exe` exits 0 |
| `W17` | `tools/check_magic_numbers.ps1` exists and runs (reports the unclassified count) |
| `W18` | `tools/check_magic_numbers.ps1 -Strict` exits 0 (i.e. **0** unclassified literals) |
| `W19` | `verify_dq_algebra.m` runs clean and `tools/check_xref.ps1` reports 0 missing |
| `W20`–`W24` | `test_integrator_rkf78/symplectic/order`, `test_precision_scaling`, `test_constants_swap` exit 0 |

Each W-gate reads **`not implemented`** until its Work-Queue item lands; that is the gate flipping to
`PASS` and the remaining-count decrementing. The next task is always the lowest-numbered `W` that is not
yet `PASS`.

Gates are driven only by live artifacts (a build result, an exe exit code, a verifier's output, a grep of
the actual source, or a checker script's exit) — never by a status word written into a doc. That keeps the
verdict honest: a gate can only pass by the thing it checks actually being true.
