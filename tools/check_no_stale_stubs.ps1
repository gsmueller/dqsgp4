<#
.SYNOPSIS
    Phase-A1 truth guard for the DQSGP4 Completion Roadmap.

.DESCRIPTION
    Fails (exit 1) if any of a set of specific stale doc comments reappears under
    src/. Each phrase below described code as unimplemented / forthcoming when it
    is in fact implemented and wired; all were corrected in Phase A1. This guard
    prevents their regression. See design/DQSGP4_COMPLETION_ROADMAP.md (Phase A).

    Scope is deliberately a fixed allow-list of KNOWN lies, not a blanket ban on
    "not yet" — honest in-progress markers (e.g. the tesseral/OMM gaps tracked as
    later roadmap items) are legitimate and must not trip this gate.
#>
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$src  = Join-Path $repo 'src'

# Literal phrases (SimpleMatch) that must NOT appear in src/.
$forbidden = @(
    'GMST polynomial. (Placeholder',          # compute_gmst IS implemented and wired
    'REQ-PR-2, forthcoming',                  # summation IS in Propagator::compute_acceleration
    '009_brouwer_rates.md',                   # dangling derivation-doc pointer (file never existed)
    'Halley Kepler solver, Aoki 1982 GMST'    # standard_sgp4() installs Newton, not Halley
)

$files = Get-ChildItem -Path $src -Recurse -Include *.h, *.hpp, *.cpp -File
$hits = @()
foreach ($phrase in $forbidden) {
    $m = $files | Select-String -SimpleMatch -Pattern $phrase
    if ($m) { foreach ($x in $m) { $hits += ("{0}:{1}" -f $x.Path, $x.LineNumber) } }
}

if ($hits.Count -gt 0) {
    Write-Output '[FAIL] check_no_stale_stubs: stale doc-lie phrase(s) present in src/:'
    $hits | ForEach-Object { Write-Output "   $_" }
    exit 1
}
Write-Output "[PASS] check_no_stale_stubs: none of $($forbidden.Count) guarded phrases present in src/"
exit 0
