<#
.SYNOPSIS
    Phase-A2 dead-scaffolding guard for the DQSGP4 Completion Roadmap.

.DESCRIPTION
    Fails (exit 1) if any orphaned / dead / phantom scaffolding removed in Phase
    A2 reappears under src/, or if the Kepler-solver de-duplication regresses
    (both model factories must delegate to the single-source orbit:: solvers
    rather than inlining their own copies). See design/DQSGP4_COMPLETION_ROADMAP.md
    (Phase A) and design/DQSGP4_ISSUE_REGISTER.md (item A2).
#>
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$src  = Join-Path $repo 'src'
$fail = @()

# 1. The orphaned KaulaTable module must be gone.
if (Test-Path (Join-Path $src 'sgp4\precomputed.h')) {
    $fail += 'src/sgp4/precomputed.h still exists (orphaned KaulaTable module)'
}

$files = Get-ChildItem -Path $src -Recurse -Include *.h, *.hpp, *.cpp -File

# 2. Dead / phantom symbols must be absent from src/.
$forbidden = @(
    'PrecomputedModel',             # phantom doc type (never existed)
    'PrecomputedSatellite',         # phantom doc type (never existed)
    'struct NearSpaceCoefficients', # unused struct
    'struct DeepSpaceCoefficients', # unused struct
    'struct PropagationState',      # unused struct
    'DragCoefficientsFn',           # broken/throwing injection slot (removed)
    '009_sgp4_modified_kepler.md'   # dangling derivation-doc pointer (file never existed)
)
foreach ($phrase in $forbidden) {
    $m = $files | Select-String -SimpleMatch -Pattern $phrase
    if ($m) { foreach ($x in $m) { $fail += ("{0}:{1}: '{2}'" -f $x.Path, $x.LineNumber, $phrase) } }
}

# 3. modified_kepler.h must not advertise the never-written Householder solver.
$mk = Join-Path $src 'orbit\modified_kepler.h'
if (Test-Path $mk) {
    foreach ($phrase in @('Householder', 'Three solver')) {
        if (Select-String -Path $mk -SimpleMatch -Pattern $phrase) {
            $fail += ("orbit/modified_kepler.h still advertises '{0}'" -f $phrase)
        }
    }
}

# 4. The Kepler de-duplication must hold: both factories delegate to the
#    single-source orbit:: solvers (A2-l).
$mf = Join-Path $src 'sgp4\model_functions.h'
$ms = Join-Path $src 'sgp4\model_selector.h'
if (-not (Select-String -Path $mf -SimpleMatch -Pattern 'orbit::solve_kepler_newton')) {
    $fail += 'model_functions.h no longer delegates to orbit::solve_kepler_newton'
}
if (-not (Select-String -Path $ms -SimpleMatch -Pattern 'orbit::solve_kepler_halley')) {
    $fail += 'model_selector.h no longer delegates to orbit::solve_kepler_halley'
}

if ($fail.Count -gt 0) {
    Write-Output '[FAIL] check_no_dead_code: dead scaffolding present or Kepler de-dup regressed:'
    $fail | ForEach-Object { Write-Output "   $_" }
    exit 1
}
Write-Output '[PASS] check_no_dead_code: no dead/orphaned scaffolding; Kepler solver single-sourced'
exit 0
