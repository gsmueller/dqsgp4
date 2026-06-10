<#
.SYNOPSIS
    W19 (companion to verify_dq_algebra.m): validate the code-to-documentation
    cross-reference is complete — every constant row is dispositioned.

.DESCRIPTION
    design/code_to_documentation_xref.md maps every reference-SGP4 constant and
    formula to its source document (or explicitly dispositions it as
    SGP4-analytical-theory-specific). This checker re-validates that mapping:
    every markdown TABLE data row that names a constant / code item must carry
    at least one non-empty disposition cell (a value, a document, a model, or a
    note). A named row with no disposition at all is "missing".

    Exit 0 iff zero missing (the W19 DoD's "check_xref reports 0 missing").
#>

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$doc  = Join-Path $repo 'design\code_to_documentation_xref.md'

if (-not (Test-Path $doc)) {
    Write-Host "  [FAIL] code_to_documentation_xref.md missing" -ForegroundColor Red
    exit 1
}

# Column-0 header labels that start a table header row (skip those, and the
# |---|---| separator rows).
$headerFirst = '^(Code|Constant|Formula|Threshold|Parameter|Item|Code Variable|Code \(both)'

$rows = 0
$missing = New-Object System.Collections.Generic.List[string]

foreach ($line in Get-Content $doc) {
    $t = $line.Trim()
    if (-not $t.StartsWith('|')) { continue }            # not a table row
    if ($t -match '^\|[\s:\-|]+\|?\s*$') { continue }     # |---|---| separator
    # inner cells (strip the outer empty cells produced by leading/trailing |)
    $cells = ($t -split '\|') | ForEach-Object { $_.Trim() }
    $cells = @($cells | Select-Object -Skip 1)            # drop leading ''
    if ($cells.Count -ge 1 -and $cells[-1] -eq '') { $cells = @($cells[0..($cells.Count-2)]) }
    if ($cells.Count -lt 2) { continue }
    $name = $cells[0]
    if ([string]::IsNullOrWhiteSpace($name)) { continue } # blank name: not a constant row
    if ($name -match $headerFirst) { continue }           # header row
    $rows++
    $disposed = @($cells[1..($cells.Count - 1)] | Where-Object { $_ -ne '' -and $_ -ne '—' })
    if ($disposed.Count -eq 0) { $missing.Add($name) }
}

Write-Host "=== check_xref: design/code_to_documentation_xref.md ===" -ForegroundColor Cyan
Write-Host ("constant rows: {0}; missing disposition: {1}" -f $rows, $missing.Count)
foreach ($m in $missing) { Write-Host "  [MISSING] $m" -ForegroundColor Yellow }

if ($missing.Count -eq 0) {
    Write-Host "xref complete: 0 missing" -ForegroundColor Cyan
    exit 0
}
exit 1
