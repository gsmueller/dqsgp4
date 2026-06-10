<#
.SYNOPSIS
    CR1-b constant-honesty scan: no decimal-truncation mis-filed as a
    measurement sigma in the SGP4 / perturbation / gravity MODEL files.

.DESCRIPTION
    CR1 (feedback_constants_generative_or_bounded) requires a finite-digit book
    value to record its decimal-truncation HONESTLY -- as a precision bound via
    from_truncated_decimal (or, for a model coefficient in the deep-space
    calling-card path, as model-fidelity accuracy via model_coefficient). It must
    NOT be smuggled in as a `measured(value, sigma)` MEASUREMENT uncertainty when
    `sigma` is really just the place value of the value's last written digit
    (the "decimal ULP"): that mislabels a representation/model limit as physical
    measurement noise.

    This scan flags every `measured("V", "S")` in the model files where S equals
    the decimal ULP of V (to within 1%). Genuine measurement sigmas -- EGM2008
    formal errors (1.3e-11, ...), the SR3 resonance sigmas (one order above the
    ULP), a published GM uncertainty (8e5) -- are NOT the ULP and pass cleanly.
    The geodesy layer is intentionally excluded: its measured() constants carry
    real published sigmas (e.g. IERS a_earth +-0.1 m) that can coincide with the
    last-digit place yet are genuine.

    Exit 0 iff no mis-filed sigma remains. Gate CR1B in run_acceptance.ps1.
#>
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot

# Model files whose measured() sigmas must be genuine (not decimal ULPs) AND
# whose finite-digit model coefficients must use model_coefficient (not
# from_truncated_decimal) so computational precision tightens with T.
$files = @(
    'src\sgp4\model_selector.h',
    'src\sgp4\deep_space.h',
    'src\perturbation\resonance.h',
    'src\constants\zonal_harmonics.h',
    'src\constants\tesseral_harmonics.h',
    'src\forces\drag.h',
    'src\forces\gravity_zonal.h',
    'src\astronomy\obliquity.h',
    'src\astronomy\earth_orbit.h',
    'src\astronomy\sidereal_time.h'
)

# Decimal ULP of a written decimal string: 10^(explicit exponent - fractional
# digits) -- the place value of the least-significant written digit. Mirrors
# TrackedValue::decimal_ulp in src/math/tracked_value.h.
function Get-DecimalUlp([string]$s) {
    $frac = 0; $exp = 0; $expSign = 1; $afterDot = $false; $inExp = $false
    foreach ($c in $s.ToCharArray()) {
        if ($c -eq '.') { $afterDot = $true }
        elseif ($c -eq 'e' -or $c -eq 'E') { $inExp = $true }
        elseif ($c -eq '-') { if ($inExp) { $expSign = -1 } }
        elseif ($c -ge '0' -and $c -le '9') {
            if ($inExp) { $exp = $exp * 10 + ([int][string]$c) }
            elseif ($afterDot) { $frac++ }
        }
    }
    $place = $expSign * $exp - $frac
    return [math]::Pow(10, $place)
}

$rx = [regex]'(?:::|\.)measured\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)'
$violations = @()

foreach ($rel in $files) {
    $path = Join-Path $repo $rel
    if (-not (Test-Path $path)) { continue }
    $lineNo = 0
    foreach ($line in (Get-Content $path)) {
        $lineNo++
        foreach ($m in $rx.Matches($line)) {
            $v = $m.Groups[1].Value
            $s = $m.Groups[2].Value
            $ulp = Get-DecimalUlp $v
            $sigma = [double]$s
            if ($ulp -gt 0 -and [math]::Abs($sigma - $ulp) / $ulp -lt 1e-6) {
                $violations += "${rel}:${lineNo}  measured(`"$v`", `"$s`") -- sigma == decimal ULP ($ulp); use from_truncated_decimal / model_coefficient"
            }
        }
    }
}

# Precision-purity check (2026-06-05 panel ruling). A finite-digit MODEL
# coefficient must NOT use from_truncated_decimal in the model files: that books
# its written-digit floor into PRECISION, which makes the computational precision
# T-INDEPENDENT and defeats the precision-tightens-with-T calling card. Model
# coefficients use model_coefficient (digit floor -> accuracy, precision scales);
# the only legitimate from_truncated_decimal here is a non-coefficient numerical
# guard (the deep_space.h 1.5e-12 retrograde-singularity regularizer), allowlisted.
$ftd = [regex]'from_truncated_decimal\(\s*"([^"]+)"'
$allow = @('1.5e-12')   # numerical guards, not model coefficients
foreach ($rel in $files) {
    $path = Join-Path $repo $rel
    if (-not (Test-Path $path)) { continue }
    $lineNo = 0
    foreach ($line in (Get-Content $path)) {
        $lineNo++
        foreach ($m in $ftd.Matches($line)) {
            $v = $m.Groups[1].Value
            if ($allow -notcontains $v) {
                $violations += "${rel}:${lineNo}  from_truncated_decimal(`"$v`") in a model file -- a finite-digit model coefficient must use model_coefficient (digit floor -> accuracy) so precision tightens with T"
            }
        }
    }
}

# defined() CONVENTION-ALLOWLIST honesty (user directive 2026-06-06: "all
# constants need accuracy and precision tracked; each J_k is the result of a
# series truncation" -- Phase 8, the strongest form). defined() may encode ONLY
# a value that is exact-by-CONVENTION (precision scales with T, accuracy 0 is
# CORRECT). Anything physical / measured / fit (every J_k, GM, eccentricity,
# obliquity, SR3 coefficient) must carry accuracy via measured()/model_coefficient;
# a transcendental (pi/4, ...) must be GENERATED so it scales with T, never
# stamped as a truncated decimal. So instead of guessing "is this physical?" by
# value range, we invert it: maintain the COMPLETE registry of accepted
# conventions and flag ANY defined() value not in it across the ENTIRE src/ tree.
# Adding a new defined() value REQUIRES registering it here with a justification
# -- that is the enforcement: every defined() is a consciously-blessed convention.
#
# Deliberate, documented exclusions (NOT flagged, by design): this scan touches
# ONLY defined(). The SR3 sgp4_standard() frozen 1970s constants (solar_system.h)
# and the IERS/geodesy presets use measured() with an adoption-bound sigma (= the
# last-digit place for a rounded adopted astronomical quantity) -- a genuine
# measurement, per the panel ruling, that must reproduce SR3 bit-for-bit (OR1);
# model_coefficient encodings carry their digit floor as accuracy. All honest.
$conventions = @(
    @{ v = '6378137.0';       why = 'WGS84/GRS80 datum semi-major axis a [m]' },
    @{ v = '6378135.0';       why = 'WGS72 datum semi-major axis a [m]' },
    @{ v = '6378.137';        why = 'WGS84/GRS80 datum a [km]' },
    @{ v = '6378.135';        why = 'WGS72 datum a [km]' },
    @{ v = '298.257223563';   why = 'WGS84 datum reciprocal flattening 1/f' },
    @{ v = '7.2921151467e-5'; why = 'WGS84/WGS72 datum Earth rotation rate omega [rad/s]' },
    @{ v = '7.292115e-5';     why = 'nominal/IERS Earth rotation rate omega [rad/s]' },
    @{ v = '2451545.0';       why = 'J2000.0 epoch Julian Date [d]' },
    @{ v = '36525.0';         why = 'Julian century in days (365.25*100, exact)' },
    @{ v = '86400.0';         why = 'SI seconds per day (exact)' },
    @{ v = '149597870700';    why = 'IAU 2012 astronomical unit AU [m] (exact by definition, Res. B2)' },
    @{ v = '299792458';       why = 'SI defining constant: speed of light c [m/s] (exact)' },
    @{ v = '3.828e26';        why = 'IAU 2015 Res. B3 nominal solar luminosity L_sun [W] (exact by convention)' },
    @{ v = '0.7790572732640'; why = 'IAU 2000 Res. B1.8 ERA defining relation: ERA(J2000) offset [rev] (exact by convention)' },
    @{ v = '1.00273781191135448'; why = 'IAU 2000 Res. B1.8 ERA defining relation: ERA rate [rev/UT1-day] (exact by convention)' }
)
$convVals = $conventions | ForEach-Object {
    $d = 0.0
    [void][double]::TryParse($_.v, [System.Globalization.NumberStyles]::Float,
                             [System.Globalization.CultureInfo]::InvariantCulture, [ref]$d)
    $d
}
$srcRoot = Join-Path $repo 'src'
$allSrc = Get-ChildItem -Path $srcRoot -Recurse -Include *.h, *.cpp |
    ForEach-Object { $_.FullName.Substring($repo.Length + 1) }
$defRx = [regex]'(?:::|\.)defined\(\s*"([^"]+)"\s*\)'
foreach ($rel in $allSrc) {
    $path = Join-Path $repo $rel
    if (-not (Test-Path $path)) { continue }
    $lineNo = 0
    foreach ($line in (Get-Content $path)) {
        $lineNo++
        foreach ($m in $defRx.Matches($line)) {
            $v = $m.Groups[1].Value
            $d = 0.0
            if (-not [double]::TryParse($v, [System.Globalization.NumberStyles]::Float,
                                        [System.Globalization.CultureInfo]::InvariantCulture, [ref]$d)) { continue }
            $ok = $false
            foreach ($cv in $convVals) {
                if ($cv -eq 0.0) { if ($d -eq 0.0) { $ok = $true; break } }
                elseif ([math]::Abs($d - $cv) -le [math]::Abs($cv) * 1e-10) { $ok = $true; break }
            }
            if (-not $ok) {
                $violations += "${rel}:${lineNo}  defined(`"$v`") is not a recognized exact-by-convention constant -- if physical/measured/fit use measured()/model_coefficient(); if transcendental GENERATE it (scales with T); if a genuinely NEW convention, register it in the conventions allowlist (tools/check_constant_honesty.ps1) with justification"
            }
        }
    }
}

if ($violations.Count -gt 0) {
    Write-Host "constant-honesty FAIL:" -ForegroundColor Yellow
    $violations | ForEach-Object { Write-Host "  $_" }
    exit 1
}
Write-Host "constant-honesty PASS: no mis-filed sigma, no precision-flooring model coefficient." -ForegroundColor Green
exit 0
