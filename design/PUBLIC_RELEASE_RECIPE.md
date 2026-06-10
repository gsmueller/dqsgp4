# Public release recipe (github.com/gsmueller/dqsgp4)

How the public repository is assembled from this (private) tree. The public repo has its
own single-commit-style history (authored `Graham Mueller <remstadt@gmail.com>`, no
tooling attribution lines, per the owner's instruction 2026-06-10); it is REPLACED, not
merged, on each release. Working area: `build/public/dqsgp4` (gitignored here).

## 1. Allowlist copy

Copy ONLY: `src/ tests/ examples/ tools/ docs/ design/ sgp4lib/ sgp4.sln vcpkg.json
README.md BUILDING.md` plus exactly three reference files (the only ones any gate reads
at runtime): `sgp4_references/aholinch_sgp4/{data/SGP4-VER.TLE, data/tcppver.out,
LICENSE}`. Never copy anything else from `sgp4_references/` (scans, third-party code
trees), the root PDFs/docx/xlsm, `archive/`, `third_party/`, `vcpkg_installed/`
(restored by the `vcpkg.json` manifest), or any build output.

## 2. The NRLMSISE-00 excision (legal: the Vallado/CelesTrak source tree is AGPL-3.0;
##    the verbatim port is a derivative and cannot ship under PolyForm-NC)

Delete in the copy: `src/atmosphere/nrlmsise00_port.h`, `src/atmosphere/nrlmsise00.h`,
`tests/test_nrlmsise00/`, `tools/gen_msis_port.py`, `tools/gen_msis_oracle.cpp`,
`design/derivations/nrlmsise00.md`. Then patch the copy for coherence:

- `sgp4.sln`: remove the `test_nrlmsise00` Project block, its four config lines, and
  its NestedProjects line (GUID search-and-delete).
- `tools/run_acceptance.ps1`: remove the MSIS1 ExeGate and its comment; reword the EX2
  comment's space-weather mention.
- `examples/quickstart.cpp`: remove the space-weather section and its include;
  renumber the following sections; update the header list and any gate-name list.
- `tools/gen_docs.py`: remove the guide's space-weather section, the MSIS symbols from
  the guide `require(...)` list, the atmosphere entry in `MODULE_USAGE`, the MSIS rows
  in the verification table and references, and `MSIS1` from the `gates_for` regex;
  reword the atmosphere module blurb.
- `tools/gen_diagrams.py`: third atmosphere-ladder rung becomes the generic
  space-weather-seam box; remove MSIS gate mentions from labels.
- `src/dqsgp4.h` and `src/forces/drag.h`: reword the NRLMSISE pointers to the generic
  DensityModel-seam language. (`nrlmsise00_density_model_stub` keeps its NAME — it is
  original Lane-fallback code with test callers; only shipped CODE was the issue.)
- Mentions of the model NAME (e.g. in measured comparisons) are fine; claims that the
  repo CONTAINS the port are not.

## 3. Public-only additions

- `LICENSE.md` = two `Required Notice:` lines (copyright Graham Mueller; commercial
  contact) + the canonical PolyForm Noncommercial 1.0.0 text
  (https://polyformproject.org/licenses/noncommercial/1.0.0.txt), verbatim.
- `THIRD-PARTY-NOTICES.md` (provenance: redistributed files, embedded data values with
  sources, algorithms-from-literature, dev-time oracles NOT redistributed, vcpkg deps,
  names). Update if data sources change.
- `.gitignore` (build outputs, `vcpkg_installed/`, `__pycache__/`).
- `.github/workflows/pages.yml` — deploys `docs/` to GitHub Pages on push
  (upload-pages-artifact + deploy-pages; the repo's Pages build_type is "workflow" —
  the legacy builder wedged).
- `README.md` additions: a `**Documentation:** <https://gsmueller.github.io/dqsgp4/>`
  pointer and the License section (PolyForm-NC, not OSI open source, commercial
  contact, third-party notices link, no-warranty/not-for-operations note).
- `BUILDING.md`: the vcpkg paragraph says manifest mode (`vcpkg install`), not
  "vendored in the tree".
- `tools/gen_docs.py`: set the published-site link base — either the `WEB_BASE`
  constant (this tree's generator has the hook) or the post-process that rewrites
  `href=…../design/derivations/` to the repository blob URL — then REGENERATE docs in
  the copy so DOC1 stays green there.
- `tests/test_sgp4/verify_deep_space_{drag,error_codes}.py`: this tree already
  contains the pip-`sgp4` fallback and graceful skip; verify they skip cleanly in the
  copy (no `sgp4_references/python_sgp4_rhodes`).

## 4. Validate the copy as its own repository

Junction `vcpkg_installed` from this tree (ignored by the copy's .gitignore), build
`sgp4.sln` from scratch, run the FULL acceptance suite in the copy. Expected: every
gate green (the suite is one smaller than this tree's, MSIS1 removed). Also run
`python tools/gen_docs.py --check` in the copy.

## 5. Publish

Fresh single-commit history in the copy (`git init -b main`, local identity
`Graham Mueller <remstadt@gmail.com>`, one release commit, NO tooling attribution
trailer), `git config http.sslBackend schannel` (this machine's TLS interception
breaks git's OpenSSL bundle), remote `https://github.com/gsmueller/dqsgp4.git`,
force-push `main`. The Pages workflow redeploys the docs automatically; verify the
site serves and spot-check pages. The repo was created and is administered via the
GitHub REST API with the GCM-stored PAT (verify the token's `/user` login is
`gsmueller` before any use; never print it).
