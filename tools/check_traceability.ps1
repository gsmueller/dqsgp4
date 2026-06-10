<#
.SYNOPSIS
    Acceptance gate W10: traceability has zero orphans.

.DESCRIPTION
    Walks the charter -> spec -> audit -> test traceability graph and fails if
    any requirement or audit identifier fails to TRACE BACK to a charter root
    (an OBJ-N or CON-N), or if any test fails to cite an audit.

    charter/objectives.md is explicit that the relation is transitive — "every
    requirement ... must trace back to an objective here." So an identifier is
    anchored if it cites a charter root directly OR cites another requirement /
    audit that is itself anchored. This admits the real, valid graph:

      REQ-IN-2  -> REQ-EF-15 -> OBJ/CON        (a spec refining another spec)
      AUD-MC-1  -> REQ-DQ-1  -> OBJ-1          (an audit verifying a spec)
      AUD-CC-1  -> CON-9, CON-13               (a constraint-conformance audit)
      AUD-DOC-2 -> AUD-CC-3  -> CON/OBJ        (an audit verifying an audit)

    An ORPHAN is a requirement or audit that cannot reach any OBJ/CON by these
    edges (it cites nothing upstream, or only participates in a free-floating
    cycle). A DANGLING citation is a concrete OBJ/CON/REQ/AUD reference whose
    target is not defined anywhere.

    Citation forms recognized in a node's own section (its heading through the
    line before the next identifier heading):
      - concrete ids: OBJ-3, CON-5, REQ-EF-15, AUD-CC-3 (any "Verifies:"/
        "Verifies"/inline form);
      - wildcard charter anchor: "every OBJ" / "all objectives" / "every CON";
      - wildcard family: "every REQ-EF-N" (edges to every defined REQ-EF-*).

    Heading levels ## .. #### are recognized; the audit layer is the five
    governance docs (first line "# Audit: ...").

    Exit 0 iff zero orphans and zero dangling citations (acceptance gate W10).
#>

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot

$fail = 0
$dangling = New-Object System.Collections.Generic.List[string]
function Orphan([string]$msg) { $script:fail++; Write-Host "  [ORPHAN] $msg" -ForegroundColor Red }

$charterFiles = Get-ChildItem (Join-Path $repo 'design\charter') -Filter *.md -File
$specFiles    = Get-ChildItem (Join-Path $repo 'design\specifications') -Filter *.md -File
$auditFiles   = @(Get-ChildItem (Join-Path $repo 'design\audit') -Filter *.md -File |
    Where-Object { (Get-Content $_.FullName -TotalCount 1) -match '^# Audit: ' })

# Carve a file into identifier sections: id -> section text (heading .. next id).
function Get-Sections($files, [string]$idRegex) {
    $map = @{}
    foreach ($f in $files) {
        $curId = $null; $buf = $null
        foreach ($line in Get-Content $f.FullName) {
            if ($line -match "^#{2,4}\s+($idRegex)\b") {
                if ($curId) { $map[$curId] = $buf.ToString() }
                $curId = $Matches[1]; $buf = [System.Text.StringBuilder]::new()
            }
            if ($curId) { [void]$buf.AppendLine($line) }
        }
        if ($curId) { $map[$curId] = $buf.ToString() }
    }
    return $map
}

# Charter roots.
$objCon = New-Object System.Collections.Generic.HashSet[string]
foreach ($f in $charterFiles) {
    foreach ($line in Get-Content $f.FullName) {
        if ($line -match '^#{2,4}\s+((?:OBJ|CON|GOAL)-\d+)\b') { [void]$objCon.Add($Matches[1]) }
    }
}

$reqSections = Get-Sections $specFiles  'REQ-[A-Z]+-\d+'
$audSections = Get-Sections $auditFiles 'AUD-[A-Z]+-\d+'
$reqIds = New-Object System.Collections.Generic.HashSet[string]; foreach ($k in $reqSections.Keys) { [void]$reqIds.Add($k) }
$audIds = New-Object System.Collections.Generic.HashSet[string]; foreach ($k in $audSections.Keys) { [void]$audIds.Add($k) }
$allDefined = New-Object System.Collections.Generic.HashSet[string]
foreach ($s in $objCon) { [void]$allDefined.Add($s) }
foreach ($s in $reqIds) { [void]$allDefined.Add($s) }
foreach ($s in $audIds) { [void]$allDefined.Add($s) }

$concreteRe = '\b((?:REQ|AUD)-[A-Z]+-\d+|(?:OBJ|CON|GOAL)-\d+)\b'
$parents = @{}                                              # id -> string[] of cited parents
$seed    = New-Object System.Collections.Generic.HashSet[string]   # ids anchored by a charter wildcard

function Process-Node([string]$id, [string]$body) {
    $ps = New-Object System.Collections.Generic.List[string]
    foreach ($m in [regex]::Matches($body, $concreteRe)) {
        $tok = $m.Groups[1].Value
        if ($tok -eq $id) { continue }                     # self-reference is not an edge
        if (-not $allDefined.Contains($tok)) { $script:dangling.Add("$id cites undefined $tok"); continue }
        [void]$ps.Add($tok)
    }
    if ($body -match '(?i)\bevery\s+OBJ\b' -or $body -match '(?i)\ball\s+objectives\b' -or
        $body -match '(?i)\bevery\s+CON\b' -or $body -match '(?i)\ball\s+constraints\b') {
        [void]$script:seed.Add($id)
    }
    foreach ($m in [regex]::Matches($body, '(?i)\bevery\s+(REQ|AUD)-([A-Z]+)-N\b')) {
        $fam  = ($m.Groups[1].Value + '-' + $m.Groups[2].Value).ToUpper()
        $pool = if ($m.Groups[1].Value.ToUpper() -eq 'REQ') { $reqIds } else { $audIds }
        foreach ($d in $pool) { if ($d -like "$fam-*") { [void]$ps.Add($d) } }
    }
    $script:parents[$id] = $ps.ToArray()
}

foreach ($id in $reqSections.Keys) { Process-Node $id $reqSections[$id] }
foreach ($id in $audSections.Keys) { Process-Node $id $audSections[$id] }

# Transitive anchoring: roots are OBJ/CON plus wildcard-seeded nodes.
$anchored = New-Object System.Collections.Generic.HashSet[string]
foreach ($x in $objCon) { [void]$anchored.Add($x) }
foreach ($x in $seed)   { [void]$anchored.Add($x) }
$changed = $true
while ($changed) {
    $changed = $false
    foreach ($id in @($parents.Keys)) {
        if ($anchored.Contains($id)) { continue }
        foreach ($p in $parents[$id]) {
            if ($anchored.Contains($p)) { [void]$anchored.Add($id); $changed = $true; break }
        }
    }
}

Write-Host "=== check_traceability: charter -> spec -> audit -> test ===" -ForegroundColor Cyan

foreach ($id in ($reqIds + $audIds)) {
    if (-not $anchored.Contains($id)) { Orphan "$id does not trace back to any OBJ/CON" }
}
foreach ($d in $dangling) { Orphan $d }

# test -> AUD: every governance test cites at least one defined audit id.
$govTests = @(Get-ChildItem (Join-Path $repo 'tests') -Filter main.cpp -Recurse -File |
    Where-Object { Select-String -Path $_.FullName -Pattern 'AUD-[A-Z]+-\d+' -Quiet })
foreach ($t in $govTests) {
    $dir = Split-Path -Leaf (Split-Path -Parent $t.FullName)
    $cited = @([regex]::Matches((Get-Content -Raw $t.FullName), 'AUD-[A-Z]+-\d+') |
        ForEach-Object { $_.Value } | Select-Object -Unique)
    if ($cited.Count -eq 0) {
        Orphan "test $dir cites no AUD identifier"
    } else {
        foreach ($c in $cited) { if (-not $audIds.Contains($c)) { Orphan "test $dir cites undefined $c" } }
    }
}

Write-Host ""
Write-Host ("graph: {0} OBJ/CON/GOAL roots, {1} REQ, {2} AUD, {3} governance tests" -f `
    $objCon.Count, $reqIds.Count, $audIds.Count, $govTests.Count) -ForegroundColor Gray
if ($fail -eq 0) {
    Write-Host "traceability: 0 orphans" -ForegroundColor Cyan
    exit 0
}
Write-Host "traceability: $fail orphan(s) / dangling citation(s)" -ForegroundColor Red
exit 1
