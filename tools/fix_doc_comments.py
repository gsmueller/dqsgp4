#!/usr/bin/env python
"""Q4 helper: promote plain `//` comment runs to `///` doc-comments EXACTLY
where the gen_docs parser reports a public entity with no attached doc and a
contiguous full-line `//` run sits directly above it. Conversion is marker-only
(the text is untouched); section banners (`// ---`) are never converted.

Usage: python tools/fix_doc_comments.py <module> [--dry]
Prints every converted run; re-run gen_docs to see the coverage change.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode = True
import gen_docs as g

PLAIN = re.compile(r"^(\s*)//($| [^/])")
BANNER = re.compile(r"^\s*// ---")


def runs_for(path, rel):
    p = g.parse_header(path, rel)
    targets = []
    for e in g.dedupe_prototypes(p.entities):
        if not e["doc"].strip():
            targets.append(e["line"])
        for m in e["members"]:
            if m["kind"] == "function" and not m["doc"].strip():
                targets.append(m["line"])
    return sorted(set(targets))


def convert(path, rel, dry):
    with open(path, encoding="utf-8", errors="replace") as f:
        head = f.read(2500)
    if "AUTO-GENERATED" in head or "DO NOT EDIT" in head.upper():
        return []  # generated files are fixed via their generators, never here
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines(keepends=True)
    changed = []
    for decl_line in runs_for(path, rel):
        # The entity's recorded line is 1-based; a doc run would END at
        # decl_line-1 (0-based index decl_line-2). A template prefix may sit
        # between: walk up over template/blank-free lines that are code.
        j = decl_line - 2
        while j >= 0 and re.match(r"\s*template\s*<", lines[j]):
            j -= 1
        if j < 0 or not PLAIN.match(lines[j]) or BANNER.match(lines[j]):
            continue
        top = j
        while top - 1 >= 0 and PLAIN.match(lines[top - 1]) and not BANNER.match(lines[top - 1]):
            top -= 1
        for k in range(top, j + 1):
            lines[k] = re.sub(r"^(\s*)//", r"\1///", lines[k], count=1)
        changed.append((top + 1, j + 1, lines[top].strip()[:70]))
    if changed and not dry:
        with open(path, "w", encoding="utf-8", newline="") as f:
            f.writelines(lines)
    return changed


def main():
    mod = sys.argv[1]
    dry = "--dry" in sys.argv
    moddir = os.path.join(g.SRC, mod)
    total = 0
    for fn in sorted(os.listdir(moddir)):
        if not fn.endswith(".h"):
            continue
        ch = convert(os.path.join(moddir, fn), f"src/{mod}/{fn}", dry)
        for a, b, first in ch:
            print(f"  {fn}:{a}-{b}  {first}")
        total += len(ch)
    print(f"{mod}: {total} run(s) {'(dry)' if dry else 'converted'}")


if __name__ == "__main__":
    main()
