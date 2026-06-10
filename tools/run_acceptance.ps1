<#
.SYNOPSIS
    The single source of truth for "is the DQSGP4 project complete?".

.DESCRIPTION
    Runs every acceptance gate defined in design/acceptance.md and prints a
    PASS / FAIL line per gate plus a final verdict: "COMPLETE" or
    "k gate(s) remaining". This script — not a human or model judgment — is the
    definition of done (gate G* in the roadmap recursive-forging-rain.md).

    Gate families:
      I0.*  regression invariant (build + test_sgp4 33/33 + module tests +
            the verify_*.m / *.py derivation verifiers) — must always hold.
      W*    one gate per Work-Queue item; most read "not implemented" until
            their item lands.

    Exit code 0 iff the verdict is COMPLETE, else 1 — so it doubles as a CI gate.

.PARAMETER SkipVerifiers
    Skip the slow Octave/Python derivation verifiers (I0.verify*, W1-5, W13-15,
    W19-dq). Use for fast iteration; a real completion run must NOT skip them.

.PARAMETER SkipBuild
    Skip the sln rebuild (assume build/Release is current).
#>
param([switch]$SkipVerifiers, [switch]$SkipBuild)

$ErrorActionPreference = 'Continue'
$repo    = Split-Path -Parent $PSScriptRoot
$msbuild = "C:\Program Files\Microsoft Visual Studio\18\Enterprise\MSBuild\Current\Bin\MSBuild.exe"
$relDir  = Join-Path $repo 'build\Release'
$gates   = [System.Collections.Generic.List[object]]::new()

# Failure signature for an Octave verifier's STDOUT (it is also failed if the
# process exit code is nonzero — Octave returns nonzero on any uncaught/parse
# error, so genuine runtime errors are caught by the exit code alone). We must
# NOT match the bare word "FAIL" or "error:" here: passing verifiers legitimately
# print "0 FAIL" / "N PASS, 0 FAIL" summaries and residual labels like
# "error:   3.3e-16" or "Relative error: = 2.4e-9". A genuine failure shows a
# NONZERO fail count, a bracketed [FAIL] status, "FAILED", or "!!!".
$failPattern = '!!!|\[FAIL\]|FAILED|[1-9][0-9]* +(FAIL|fail)'

function Gate($id, $label, [bool]$pass, $detail = '') {
    $gates.Add([pscustomobject]@{ Id = $id; Label = $label; Pass = $pass; Detail = $detail })
}

# Resolve an Octave CLI if one is available.
function Find-Octave {
    $known = 'C:\Program Files\GNU Octave\Octave-11.1.0\mingw64\bin\octave-cli.exe'
    if (Test-Path $known) { return $known }
    $c = Get-Command octave-cli -ErrorAction SilentlyContinue
    if ($c) { return $c.Source }
    return $null
}

Write-Host "=== DQSGP4 acceptance run ===" -ForegroundColor Cyan

# --- I0.build -------------------------------------------------------------
if (-not $SkipBuild) {
    & $msbuild (Join-Path $repo 'sgp4.sln') /p:Configuration=Release /p:Platform=x64 /m /verbosity:quiet *> $null
    Gate 'I0.build' 'sln builds (Release|x64)' ($LASTEXITCODE -eq 0) "MSBUILD_EXIT=$LASTEXITCODE"
} else {
    Gate 'I0.build' 'sln builds (skipped)' $true 'skipped'
}

# --- I0.sgp4 : the 33/33 regression --------------------------------------
$sgp4Exe = Join-Path $relDir 'test_sgp4.exe'
if (Test-Path $sgp4Exe) {
    $o = & $sgp4Exe 2>&1
    $sat = [bool]($o | Select-String 'Satellites:\s*33 pass,\s*0 fail')
    $pts = [bool]($o | Select-String 'Data points:\s*623 pass,\s*0 fail')
    Gate 'I0.sgp4' 'test_sgp4 = 33/33, 623/623' ($sat -and $pts) ''
} else { Gate 'I0.sgp4' 'test_sgp4 = 33/33, 623/623' $false 'exe missing' }

# --- I0.<module> : every module test exits 0 -----------------------------
$modules = 'test_math','test_geodesy','test_wgs84','test_astronomy','test_perturbation',
           'test_tle','test_propagator','test_dual_number','test_quaternion','test_dual_quaternion'
foreach ($m in $modules) {
    $exe = Join-Path $relDir "$m.exe"
    if (Test-Path $exe) { & $exe *> $null; Gate "I0.$m" "$m exits 0" ($LASTEXITCODE -eq 0) '' }
    else { Gate "I0.$m" "$m exits 0" $false 'exe missing' }
}

# --- I0.verify : Octave verify_*.m + Python verify_*.py ------------------
if ($SkipVerifiers) {
    Gate 'I0.verify' 'derivation verifiers (skipped)' $true 'skipped (-SkipVerifiers)'
} else {
    $octave = Find-Octave
    if (-not $octave) {
        Gate 'I0.verify' 'Octave verify_*.m' $false 'octave-cli not found'
    } else {
        $mfiles = Get-ChildItem -Path (Join-Path $repo 'design') -Recurse -Filter 'verify_*.m'
        $failed = @()
        foreach ($f in $mfiles) {
            $out = & $octave --no-gui --quiet $f.FullName 2>&1
            if (($LASTEXITCODE -ne 0) -or ($out | Select-String -CaseSensitive $failPattern)) {
                $failed += $f.Name
            }
        }
        Gate 'I0.verify' "Octave verify_*.m ($($mfiles.Count) files)" ($failed.Count -eq 0) `
             ($(if ($failed.Count) { "failed: $($failed -join ', ')" } else { 'all pass' }))
    }
    $py = Get-Command python -ErrorAction SilentlyContinue
    $pyfiles = Get-ChildItem -Path (Join-Path $repo 'tests\test_sgp4') -Filter 'verify_*.py' -ErrorAction SilentlyContinue
    if ($py -and $pyfiles) {
        $pyFailed = @()
        foreach ($f in $pyfiles) {
            Push-Location $repo; & $py $f.FullName *> $null; $code = $LASTEXITCODE; Pop-Location
            if ($code -ne 0) { $pyFailed += $f.Name }
        }
        Gate 'I0.pyverify' "Python verify_*.py ($($pyfiles.Count) files)" ($pyFailed.Count -eq 0) `
             ($(if ($pyFailed.Count) { "failed: $($pyFailed -join ', ')" } else { 'all pass' }))
    }
}

# --- W-gates : one per Work-Queue item -----------------------------------
# A verifier gate: PASS iff the .m exists and runs clean (else "not implemented").
function VerifierGate($id, $relPath) {
    $p = Join-Path $repo $relPath
    if (-not (Test-Path $p)) { Gate $id (Split-Path $relPath -Leaf) $false 'not implemented'; return }
    if ($SkipVerifiers) { Gate $id (Split-Path $relPath -Leaf) $false 'skipped (-SkipVerifiers)'; return }
    $octave = Find-Octave
    if (-not $octave) { Gate $id (Split-Path $relPath -Leaf) $false 'octave-cli not found'; return }
    $out = & $octave --no-gui --quiet $p 2>&1
    $ok = ($LASTEXITCODE -eq 0) -and -not ($out | Select-String -CaseSensitive $failPattern)
    Gate $id (Split-Path $relPath -Leaf) $ok ''
}
# An exe gate: PASS iff build/Release/<name>.exe exists and exits 0.
function ExeGate($id, $name) {
    $exe = Join-Path $relDir "$name.exe"
    if (-not (Test-Path $exe)) { Gate $id "$name exits 0" $false 'not implemented'; return }
    & $exe *> $null; Gate $id "$name exits 0" ($LASTEXITCODE -eq 0) ''
}
# A grep gate: PASS iff a file exists and matches a pattern.
function GrepGate($id, $label, $relPath, $pattern) {
    $p = Join-Path $repo $relPath
    if (-not (Test-Path $p)) { Gate $id $label $false 'file missing'; return }
    Gate $id $label ([bool](Select-String -Path $p -Pattern $pattern)) ''
}
# A script gate: PASS iff a checker .ps1 exists and exits 0.
function ScriptGate($id, $label, $relPath, [string[]]$scriptArgs = @()) {
    $p = Join-Path $repo $relPath
    if (-not (Test-Path $p)) { Gate $id $label $false 'not implemented'; return }
    & pwsh -NoProfile -File $p @scriptArgs *> $null; Gate $id $label ($LASTEXITCODE -eq 0) ''
}

# P1 fidelity (W1-W5)
VerifierGate 'W1' 'design\derivations\verify_near_short_period.m'
VerifierGate 'W2' 'design\derivations\verify_near_longperiod.m'
VerifierGate 'W3' 'design\derivations\verify_near_kepler.m'
VerifierGate 'W4' 'design\derivations\verify_near_secular.m'
VerifierGate 'W5' 'design\derivations\verify_near_osculating.m'
# Track B accuracy (W6-W8)
GrepGate 'W6' 'gravity_central writes errors.accuracy' 'src\forces\gravity_central.h' 'errors\.accuracy\s*='
GrepGate 'W7' 'gravity_zonal writes errors.accuracy'  'src\forces\gravity_zonal.h'   'errors\.accuracy\s*='
ExeGate  'W8' 'test_force_models'
# Track C traceability (W9-W12)
ScriptGate 'W9'  'index.md synced to reality'        'tools\check_index.ps1'
ScriptGate 'W10' 'traceability: 0 orphans'           'tools\check_traceability.ps1'
ExeGate    'W11' 'test_code_consistency'
ExeGate    'W12' 'test_error_framework'
# A7 dimensional audits (W13-W15)
VerifierGate 'W13' 'design\derivations\verify_dim_brouwer.m'
VerifierGate 'W14' 'design\derivations\verify_dim_kaula.m'
VerifierGate 'W15' 'design\derivations\verify_dim_recovery.m'
# Track D quality (W16-W19)
ExeGate    'W16' 'test_model_value'
ScriptGate 'W17' 'magic-number scanner runs'         'tools\check_magic_numbers.ps1'
# W18 = the magic-number scanner in strict mode (exit 0 iff zero unclassified).
ScriptGate 'W18' 'zero unclassified literals' 'tools\check_magic_numbers.ps1' @('-Strict')
VerifierGate 'W19' 'design\derivations\verify_dq_algebra.m'
ScriptGate   'W19' 'xref: 0 missing'                 'tools\check_xref.ps1'
# Track E future-proofing (W20-W24)
ExeGate 'W20' 'test_integrator_rkf78'
ExeGate 'W21' 'test_integrator_symplectic'
ExeGate 'W22' 'test_integrator_order'
# RK1: the unified ButcherTableau + rk_step driver (euler/rk4/rkf78 = one driver,
# three tableaux). Theory: design/derivations/runge_kutta_lie_group.md.
ExeGate 'RK1' 'test_butcher_tableau'
ExeGate 'W23' 'test_precision_scaling'
ExeGate 'W24' 'test_constants_swap'

# --- Frozen SGP4 self-consistency oracle (TEMPORARY — redesign regression) ----
# A golden-master far tighter than the 33/33 reference test; catches sub-tolerance
# drift in the SGP4 path while the common-library redesign (esp. F3) is underway.
ExeGate 'OR1' 'test_sgp4_regression'

# === DQSGP4 Completion Roadmap (design/DQSGP4_COMPLETION_ROADMAP.md) ======
# Tracked by design/DQSGP4_ISSUE_REGISTER.md. Gates read "not implemented"
# until their item lands, exactly like the W-series; victory = all PASS.
# Phase A — truthful baseline (R1).
ScriptGate   'A1'   'no stale doc-lie comments'    'tools\check_no_stale_stubs.ps1'
ScriptGate   'A2'   'no dead/orphaned scaffolding' 'tools\check_no_dead_code.ps1'
# Phase B — precision infrastructure (tracked transcendentals).
ExeGate      'B1'   'test_tracked_exp_log'
ExeGate      'B2'   'test_tracked_pow_cbrt'
ExeGate      'B3'   'test_tracked_transcendentals'
# Honest constant representation (generative or decimal-truncation-bounded).
ExeGate      'CR1'  'test_constant_representation'
# Phase C — earth models / Astronomical Almanac.
ExeGate      'C1'   'test_earth_constants'
ExeGate      'C2'   'test_earth_models'
# Phase D / Generative Astro Library — geopotential extension.
ExeGate      'GAL1' 'test_zonal_harmonics'
ExeGate      'D1'   'test_gravity_zonal_jn'
# Phase E — updated TLE / OMM ingestion.
ExeGate      'E1'   'test_tle_alpha5'
ExeGate      'E2'   'test_omm_kvn'
# Phase F — API parity.
ExeGate      'F1'   'test_state_conversion'
ExeGate      'F2'   'test_api_parity'
ExeGate      'F3'   'test_constants_single_source'
# Phase G — authentic vs boosted mode.
ExeGate      'G1'   'test_modes'
# R3a — DQ facade force injection + DqForceOptions presets (lunisolar/drag/SRP;
# default model bit-unchanged). Theory: design/derivations/dq_propagator_facade.md.
ExeGate      'FM1'  'test_force_presets'
# Cross-cutting — correctness & precision.
ExeGate      'AD1'  'test_adaptive_stepping'
ExeGate      'DS1'  'test_sdp4_precision'
ExeGate      'BUG1' 'test_celestial_body_fields'
# Remaining roadmap items (D2 tesseral, E3 OMM-XML, EPH, DRAG1, INJ1, H1) are
# ruled in-scope (design/DQSGP4_EXECUTION_PLAN.md) and gated here as they land.
ExeGate      'DRAG1' 'test_drag_density'
# R1 — piecewise-exponential static atmosphere (Vallado 8-4 / USSA76-CIRA, in-repo
# born-digital ATMOSEXP.DAT): node-exact + chain-verified table DensityModel.
# Theory: design/derivations/atmosphere_exponential_table.md.
ExeGate      'ATM1'  'test_atmosphere'
# R2 — solar radiation pressure: cannonball + cylindrical shadow; P_1AU GENERATED
# from L_sun/(4*pi*AU^2*c) (IAU B3 nominal + IAU 2012 AU + SI c, all defined).
# Theory: design/derivations/solar_radiation_pressure.md.
ExeGate      'SRP1'  'test_srp'
ExeGate      'EPH'   'test_ephemeris'
ExeGate      'D2'    'test_tesseral'
ExeGate      'INJ1'  'test_injection'
ExeGate      'E3'    'test_omm_xml'
ExeGate      'H1'    'test_attitude_dynamics'
# Series-based constants (user directive 2026-06-05): non-defined constants
# carry series-generated precision (tightens with T) + series-truncation
# accuracy (tightens with kept terms). Obliquity ε_A(t) is the exemplar.
ExeGate      'SC1'   'test_series_constants'
# Professional-library re-architecture (design/PROFESSIONAL_LIBRARY_PLAN.md).
# L1 — time scales & epochs (Epoch/TimeScale + views). Theory:
# design/derivations/time_scales_and_epochs.md.
ExeGate      'TIME1' 'test_time_scales'
# L2 — reference frames (rot_x/y/z + sidereal_rotation = Rz(GMST), Matrix3).
# Theory: design/derivations/reference_frames.md.
ExeGate      'FRAME1' 'test_reference_frames'
# L3 frame chain — ecliptic-of-date -> GCRS: IAU 2006 bias-precession (FW angles)
# + obliquity rotation Rx(-eps_A), wiring the obliquity generator. Element-wise vs
# erfa.pmat06/ecm06/obl06 (conformance). Theory: design/derivations/frame_chain.md.
ExeGate      'FRAME2' 'test_frame_chain'
# R4b — IAU 2000A luni-solar nutation (678-term in-repo SOFA table, planetary
# floors) + ERA/GMST06/GAST + polar motion + the full GCRS->ITRS chain, all
# erfa-arbitrated (nut06a/gmst06/gst06a/pom00/c2t06a; pom BIT-exact, chain
# mas-grade). Theory: design/derivations/nutation_itrs.md.
ExeGate      'NUT1'  'test_nutation'
# R4a — the gated ground-track example: TLE -> DQSGP4(+lunisolar) -> ITRS ->
# subpoint, umbrella-header-only, physical bounds asserted (examples cannot rot).
ExeGate      'EX1'   'example_ground_track'
# Q1 — the gated quickstart tour: both propagators, presets, adaptive stepping,
# the three-error budget (incl. the conservative DQ LTE-proxy semantics), T
# scaling — the source the docs guide extracts from.
ExeGate      'EX2'   'example_quickstart'
# Formula-layer Stage 2 — unified geopotential: monopole + zonal + tesseral in ONE
# Cunningham V/W pass (GravityField fuses the zonal column ⊕ tesseral block).
# Round-off (NOT bit) vs summed legacy forces + closed-form J2 + monopole identity.
# Theory: design/derivations/geopotential.md.
ExeGate      'GEOPOT' 'test_geopotential'
# L4 — Cartesian third-body perturbation force (Sun/Moon). Battin f(q) stable form
# + the L3 ephemeris; gated vs the naive form @bf50 [formula] and JPL DE430 body
# positions [ephemeris->accel]. Theory: design/derivations/third_body_perturbation.md.
ExeGate      'TB1'   'test_third_body'
ScriptGate   'CR1B'  'constant honesty (no over-claimed decimals)' 'tools\check_constant_honesty.ps1'
# Q1 — docs freshness: docs/ must equal what gen_docs.py emits from the current
# tree (API enumeration parsed strictly, guide snippets extracted from gated
# sources, guide symbols verified against the parsed model).
ScriptGate   'DOC1'  'generated docs fresh (regenerate-and-diff)' 'tools\check_docs_fresh.ps1'
# V1 (vision review, user decision 2026-06-10) — the DQ solver core UNPEELS from
# the retained SGP4 oracle tier: outside src/sgp4/, only the named TLE adapters
# (state_from_tle, state_conversion, orbit/state_from_elements) include sgp4/.
ScriptGate   'ARCH1' 'DQ core unpeelable (sgp4 includes contained)' 'tools\check_sgp4_containment.ps1'
# V2 — the solver's own arbitrary-precision gate (the DS1 analogue): a 600 s
# Cartesian-seeded DQ arc at double vs cpp_bin_float_50 — values agree < 1 mm;
# with the series/iteration tolerance scaled with T (1e-12 -> 1e-40) the
# position/twist PRECISION channels tighten >= 1e20 (measured ~1e26); at the
# double-grade tolerance the J2 e²-iteration floor dominates identically at
# both T (pinned: T alone is not the dial — (T, tolerance) is). The TU includes
# no sgp4/ or tle/ header: the compile itself demonstrates the unpeel.
ExeGate      'PREC1' 'test_dq_precision'

# --- report ---------------------------------------------------------------
Write-Host ""
foreach ($g in $gates) {
    $tag = if ($g.Pass) { 'PASS' } else { 'FAIL' }
    $col = if ($g.Pass) { 'Green' } else { 'Yellow' }
    $line = "  [{0}] {1,-10} {2}" -f $tag, $g.Id, $g.Label
    if ($g.Detail) { $line += "  ($($g.Detail))" }
    Write-Host $line -ForegroundColor $col
}
$passN = ($gates | Where-Object Pass).Count
$total = $gates.Count
$remaining = $total - $passN
Write-Host ""
Write-Host ("  {0}/{1} gates pass" -f $passN, $total)
if ($remaining -eq 0) {
    Write-Host "  VERDICT: COMPLETE" -ForegroundColor Green
    exit 0
} else {
    Write-Host ("  VERDICT: {0} gate(s) remaining" -f $remaining) -ForegroundColor Yellow
    exit 1
}
