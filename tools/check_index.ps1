<#
.SYNOPSIS
    Acceptance gate W9: assert design/index.md is synchronized with disk.

.DESCRIPTION
    The index's "Coverage status" table claims a file count for each
    governance spine layer (Charter -> Specifications -> Audit -> Tests).
    This script recomputes every count directly from the files on disk and
    fails if the table is stale; it additionally fails if any specification
    or audit file is not linked from the index.

    Disk-count rule, one per spine layer:
      Charter        design/charter/*.md
      Specifications design/specifications/*.md
      Audit          design/audit/*.md whose first line is "# Audit: ..."
                     (the governance audit docs; the SGP4-derivation working
                     notes that also live under audit/ lack that heading and
                     are excluded, as are the audit/ subdirectories)
      Tests          tests/*/main.cpp citing at least one AUD-<area>-N ID

    The "Module designs" row is the open narrative layer and is intentionally
    not counted.

    Exit 0 iff every check passes, so it doubles as acceptance gate W9
    (invoked by tools/run_acceptance.ps1).
#>

$ErrorActionPreference = 'Stop'
$repo  = Split-Path -Parent $PSScriptRoot
$index = Join-Path $repo 'design\index.md'

$fail = 0
function Check([string]$label, [bool]$ok, [string]$detail = '') {
    if ($ok) {
        Write-Host "  [PASS] $label" -ForegroundColor Green
    } else {
        $script:fail++
        Write-Host "  [FAIL] $label  $detail" -ForegroundColor Red
    }
}

Write-Host "=== check_index: design/index.md vs disk ===" -ForegroundColor Cyan

if (-not (Test-Path $index)) {
    Write-Host "  [FAIL] design/index.md missing" -ForegroundColor Red
    exit 1
}
$indexText = Get-Content -Raw $index

# --- recompute spine-layer counts from disk ------------------------------
$charterCount = (Get-ChildItem (Join-Path $repo 'design\charter') -Filter *.md -File).Count
$specCount    = (Get-ChildItem (Join-Path $repo 'design\specifications') -Filter *.md -File).Count

# Governance audit docs are the top-level audit/*.md whose first line is the
# "# Audit: <area>" heading. Non-recursive Get-ChildItem already skips the
# theoretical_basis_audit/ and remediation/ subdirectories.
$auditFiles = @(Get-ChildItem (Join-Path $repo 'design\audit') -Filter *.md -File |
    Where-Object { (Get-Content -Path $_.FullName -TotalCount 1) -match '^# Audit: ' })
$auditCount = $auditFiles.Count

# Governance tests are the test programs whose source cites an AUD-<area>-N ID.
$testFiles = @(Get-ChildItem (Join-Path $repo 'tests') -Filter main.cpp -Recurse -File |
    Where-Object { Select-String -Path $_.FullName -Pattern 'AUD-(CC|EF|MC|DOC|TC)-\d+' -Quiet })
$testCount = $testFiles.Count

# --- parse the coverage-table claim for one layer ------------------------
function Get-TableCount([string]$layer) {
    $m = [regex]::Match($indexText, "(?m)^\|\s*$layer\s*\|\s*(\d+)\s*/\s*(\d+)\s*\|")
    if (-not $m.Success) { return $null }
    return [pscustomobject]@{ Num = [int]$m.Groups[1].Value; Den = [int]$m.Groups[2].Value }
}

$layers = @(
    @{ Name = 'Charter';        Disk = $charterCount },
    @{ Name = 'Specifications'; Disk = $specCount },
    @{ Name = 'Audit';          Disk = $auditCount },
    @{ Name = 'Tests';          Disk = $testCount }
)

foreach ($L in $layers) {
    $claim = Get-TableCount $L.Name
    if ($null -eq $claim) {
        Check "$($L.Name): coverage row present and numeric" $false "no '| $($L.Name) | n / m |' row found"
        continue
    }
    $ok = ($claim.Num -eq $L.Disk) -and ($claim.Den -eq $L.Disk)
    Check "$($L.Name): table $($claim.Num) / $($claim.Den) equals disk count $($L.Disk)" $ok `
          "table says $($claim.Num) / $($claim.Den), disk has $($L.Disk)"
}

# --- every specification file is linked from the index -------------------
foreach ($f in Get-ChildItem (Join-Path $repo 'design\specifications') -Filter *.md -File) {
    Check "specification linked: $($f.Name)" ($indexText.Contains("specifications/$($f.Name)")) 'not linked in index'
}

# --- every governance audit file is linked from the index ----------------
foreach ($f in $auditFiles) {
    Check "audit linked: $($f.Name)" ($indexText.Contains("audit/$($f.Name)")) 'not linked in index'
}

# --- every governance test is referenced in the test plan ----------------
foreach ($f in $testFiles) {
    $dir = Split-Path -Leaf (Split-Path -Parent $f.FullName)
    Check "test referenced: $dir" ($indexText.Contains($dir)) 'not referenced in section 5'
}

Write-Host ""
if ($fail -eq 0) {
    Write-Host "index.md in sync: $charterCount charter, $specCount spec, $auditCount audit, $testCount test" -ForegroundColor Cyan
    exit 0
}
Write-Host "$fail index check(s) failed" -ForegroundColor Red
exit 1
