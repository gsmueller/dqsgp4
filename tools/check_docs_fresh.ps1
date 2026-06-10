# check_docs_fresh.ps1 — DOC1: the generated documentation set must equal what
# tools/gen_docs.py emits from the CURRENT tree (regenerate-and-diff freshness).
# Also fails on parser warnings (an unclassifiable construct in a public header)
# and on guide-marker / guide-symbol drift, since gen_docs --check enforces both.
$repo = Split-Path $PSScriptRoot -Parent
Push-Location $repo
try {
    & python (Join-Path $PSScriptRoot 'gen_docs.py') --check
    exit $LASTEXITCODE
} finally {
    Pop-Location
}
