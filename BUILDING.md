# Building & testing

Header-only C++20 library. The only compiled units are `src/tle/tle_parser.cpp` and the test `main.cpp`
files. Everything else is templates included from `src/`.

## Prerequisites

- **Visual Studio 2026** (MSVC v14.5x toolset, C++20 / `/std:c++20`).
- **vcpkg dependencies** — Boost (multiprecision + math), declared in `vcpkg.json`.
  From the repo root run `vcpkg install` (manifest mode); it materializes the expected
  headers under `vcpkg_installed/x64-windows/include/`. (Visual Studio with vcpkg
  integration does this automatically on first build.)

## Option 1 — Visual Studio (simplest)

Open `sgp4.sln`, select **Release | x64**, and Build Solution. All targets (`sgp4lib` + the `test_*`
projects) build to `build/Release/`. This always works because the IDE locates its own toolset.

## Option 2 — Command line (MSBuild)

Invoke MSBuild by **full path** — it locates the VC toolset itself, so no `vcvarsall` is needed:

```powershell
& "C:\Program Files\Microsoft Visual Studio\18\Enterprise\MSBuild\Current\Bin\MSBuild.exe" `
    "sgp4.sln" /p:Configuration=Release /p:Platform=x64 /m
```

Build a single target with `/t:` (e.g. `/t:test_sgp4`). Outputs land in `build/Release/`.

> **Gotcha — `MSBUILD_EXIT=9009`.** The legacy `build.bat` calls `vcvarsall.bat x64` then a bare
> `msbuild`. On some setups `vcvarsall` fails to put `msbuild` on `PATH` (exit 9009, "not recognized").
> The full-path invocation above sidesteps this. (`vcvarsall`/`vcvars64` *do* still set up `cl.exe`
> correctly — only the `msbuild` PATH entry was affected.)
>
> **Gotcha — LF line endings.** `build.bat` / `build_test.bat` are stored with LF endings; `cmd.exe`
> mis-parses them. If you script the build, write your temp `.bat` with **CRLF** endings.

## Option 3 — Compile a single file with `cl` (quick syntax check)

The three `src/` drivers (`main.cpp`, `test_sgp4_ver.cpp`, `test_series.cpp`) are not in the solution. To
compile-check any file:

```bat
call "C:\Program Files\Microsoft Visual Studio\18\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
cl /nologo /std:c++20 /EHsc ^
   /I vcpkg_installed\x64-windows\include /I src ^
   /c src\test_series.cpp
```

`build_test.bat` does the same to produce a standalone `test_sgp4.exe` at the repo root (separate from the
solution's `build/Release/test_sgp4.exe`).

## Running the tests

Run the executables from `build/Release/` **with the repository root as the working directory** — the
SGP4 regression reads its reference data from repo-relative paths
(`sgp4_references/aholinch_sgp4/data/SGP4-VER.TLE` and `tcppver.out`):

```powershell
build\Release\test_sgp4.exe          # headline regression: 33/33 satellites, 623/623 points
build\Release\test_math.exe          # per-module unit suites
build\Release\test_propagator.exe    # DQSGP4 LEO smoke test
# ... test_geodesy, test_wgs84, test_astronomy, test_perturbation, test_tle,
#     test_dual_number, test_quaternion, test_dual_quaternion
```

Each test exe returns 0 on success and prints a pass/fail summary.

## Derivation verifiers (optional)

The SGP4 coefficient derivations are symbolically verified in Octave/SymPy. With Octave + the symbolic
package installed, run any `design/derivations/verify_*.m` to re-check a `simplify(derived−code)=0`
identity. The deep-space coefficient checks are Python (`tests/test_sgp4/verify_deep_space_*.py`).
