# check_sgp4_containment.ps1 — ARCH1: the DQ solver core UNPEELS from the SGP4
# oracle tier. Outside src/sgp4/ itself, only the named TLE adapters may
# include a sgp4/ header (user decision 2026-06-10,
# design/DQSGP4_VISION_REVIEW.md §5: the arbitrary-precision SGP4 propagator
# is retained in full; this gate pins the layering, not a removal).
$repo = Split-Path $PSScriptRoot -Parent
$src = Join-Path $repo 'src'

# The complete sanctioned list (measured topology at gate creation):
#   - dynamics/state_from_tle.h     the TLE -> epoch State adapter (SGP4 at t=0)
#   - dynamics/state_conversion.h   the State <-> sgp4::StateVector adapters
#   - orbit/state_from_elements.h   oracle-tier pipeline piece (SGP4-internal)
$allow = @(
    'dynamics\state_from_tle.h',
    'dynamics\state_conversion.h',
    'orbit\state_from_elements.h'
)

$violations = @()
Get-ChildItem -Path $src -Recurse -Include *.h, *.cpp | ForEach-Object {
    $rel = $_.FullName.Substring($src.Length + 1)
    if ($rel -like 'sgp4\*') { return }          # the tier itself
    if ($allow -contains $rel) { return }        # the sanctioned adapters
    $hits = Select-String -Path $_.FullName -Pattern '#include\s+"(\.\./)*sgp4/'
    foreach ($h in $hits) {
        $violations += "src\$rel`:$($h.LineNumber): $($h.Line.Trim())"
    }
}

if ($violations.Count -gt 0) {
    Write-Host 'ARCH1 FAIL — sgp4/ included outside the sanctioned adapters:'
    $violations | ForEach-Object { Write-Host "  $_" }
    exit 1
}
Write-Host 'ARCH1 PASS — the DQ core unpeels (sgp4/ includes contained to the adapters).'
exit 0
