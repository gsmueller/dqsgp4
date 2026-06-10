<#
.SYNOPSIS
    Acceptance gates W17 / W18: scan the reusable library for magic numbers.

.DESCRIPTION
    Implements the zero-magic-numbers policy (design/zero_magic_numbers_policy.md,
    CON-1 / AUD-CC-15) for the REUSABLE LIBRARY — the five directories the audit
    governs: src/math, src/dynamics, src/constants, src/forces, src/integrators.

    Scope rationale. AUD-CC-15 / code_consistency.md bind the no-magic-numbers
    rule to those five directories. The SGP4 *application* code (src/sgp4,
    src/perturbation, src/orbit, src/atmosphere, src/ephemeris, src/astronomy,
    src/geodesy) reproduces the authentic published SGP4 algorithm; its numeric
    constants ARE the algorithm (governed by the 33/33 reference-fidelity gate),
    not opaque library literals, so they are deliberately out of scope here.
    src/dynamics/state_from_tle.h (the SGP4 bridge) is excluded for the same
    reason it is excluded from test_code_consistency.

    What counts. A floating-point literal (decimal point or exponent, optional
    f/l suffix) appearing in CODE — not in a // or /* */ comment, not inside a
    "..." string (defining parameters use the exact string form T("..."), which
    is therefore never flagged). Pure integer literals are not magic (AUD-CC-15
    permits small integers as exact<T>(n) arguments).

    Classification. A flagged literal is CLASSIFIED (allowed) iff its line
    carries a "magic-ok" justification comment naming why the literal is
    irreducible (e.g. the 1.0 mantissa of ldexp(1.0,n) = 2^n). Everything else
    is UNCLASSIFIED and should be rewritten to an exact form (exact<T>(n),
    T(p)/T(q), T("...")) or to a ModelValue<T> / CurveFit<T> with provenance.

.PARAMETER Strict
    Exit nonzero if any unclassified literal remains (gate W18). Without it the
    scanner always exits 0 after listing the count (gate W17).
#>

param([switch]$Strict)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$dirs = 'math', 'dynamics', 'constants', 'forces', 'integrators' |
    ForEach-Object { Join-Path $repo "src\$_" }

# A C++ floating-point literal: needs a decimal point or an exponent (pure
# integers are excluded); optional f/l suffix; not glued to an identifier/dot.
$floatRe = '(?<![\w.])(?:[0-9]*\.[0-9]+(?:[eE][-+]?[0-9]+)?|[0-9]+[eE][-+]?[0-9]+)[fFlL]?(?![\w.])'

$total = 0
$unclassified = New-Object System.Collections.Generic.List[string]
$classified = 0

foreach ($dir in $dirs) {
    if (-not (Test-Path $dir)) { continue }
    foreach ($f in Get-ChildItem $dir -Recurse -Include *.h, *.cpp -File) {
        if ($f.Name -eq 'state_from_tle.h') { continue }   # SGP4 application bridge
        $n = 0
        foreach ($line in Get-Content $f.FullName) {
            $n++
            $t = $line.TrimStart()
            if ($t -match '^(///|//|\*|/\*)') { continue }  # full-line comment
            $code = $line -replace '//.*$', ''              # drop trailing comment
            $code = $code -replace '"[^"]*"', ''            # drop string contents
            foreach ($m in [regex]::Matches($code, $floatRe)) {
                $total++
                if ($line -match 'magic-ok') {
                    $classified++
                } else {
                    $unclassified.Add(("{0}:{1}: {2}" -f $f.Name, $n, $m.Value))
                }
            }
        }
    }
}

Write-Host "=== check_magic_numbers: src/{math,dynamics,constants,forces,integrators} ===" -ForegroundColor Cyan
Write-Host ("float literals in code: {0} total, {1} classified (magic-ok), {2} unclassified" -f `
    $total, $classified, $unclassified.Count)
foreach ($u in $unclassified) { Write-Host "  [UNCLASSIFIED] $u" -ForegroundColor Yellow }

if ($Strict) {
    if ($unclassified.Count -eq 0) {
        Write-Host "STRICT: zero unclassified literals" -ForegroundColor Green
        exit 0
    }
    Write-Host "STRICT: $($unclassified.Count) unclassified literal(s) remain" -ForegroundColor Red
    exit 1
}
exit 0
