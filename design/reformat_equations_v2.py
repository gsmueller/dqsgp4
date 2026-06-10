"""
Reformat equations in Brouwer-Hori 1961 markdown - Phase 2.
Split long lines at major + or - boundaries by inserting \\\\ \\n&\\quad.
Also split long \\left[...\\right] with \\times markers.

Strategy: Find lines within $$ blocks that are > 300 chars and
don't already have \\\\ at the end. For each, find a good split point
near the midpoint at a + or - sign that's at the top level (not inside braces).
"""

import re

path = r"sgp4_references\vallado_celestrak\documentation\SGP4\Brouwer_Hori 1961\Brouwer_Hori_1961a_Atmospheric_Drag_AJ66_193.md"

with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

lines = content.split('\n')

# Count how many lines are still > 300 chars and contain LaTeX
count_before = sum(1 for l in lines if len(l.rstrip()) > 300 and '\\frac' in l)
print(f"Lines > 300 chars with LaTeX before: {count_before}")

# For now, let's just report what we'd change rather than auto-modify,
# since these structural changes need careful handling.
# The remaining 300+ char lines are mostly coefficient lines within aligned blocks.

# Let me identify the structural patterns:
in_eq = False
eq_start = -1
long_eq_lines = []
for i, line in enumerate(lines):
    s = line.rstrip()
    if '$$' in s:
        dollar_positions = [m.start() for m in re.finditer(r'\$\$', s)]
        for pos in dollar_positions:
            in_eq = not in_eq
            if in_eq:
                eq_start = i

    if in_eq and len(s) > 300 and '\\frac' in s:
        # Count brace depth to find top-level + or - signs
        depth = 0
        split_candidates = []
        for j, ch in enumerate(s):
            if ch in '{[(': depth += 1
            elif ch in '}])': depth -= 1
            elif ch in '+-' and depth == 0 and j > 50 and j < len(s) - 50:
                split_candidates.append(j)

        if split_candidates:
            # Find split point nearest to middle
            mid = len(s) // 2
            best = min(split_candidates, key=lambda x: abs(x - mid))
            long_eq_lines.append((i+1, len(s), len(split_candidates), best))

print(f"\nLines that could be split (>300c, has top-level +/- near middle):")
for ln, length, n_candidates, best_pos in long_eq_lines[:20]:
    print(f"  L{ln} ({length}c): {n_candidates} split candidates, best at char {best_pos}")

print(f"\nTotal candidates: {len(long_eq_lines)}")
print("\nThese lines are dense coefficient expansions within aligned blocks.")
print("They could be split by inserting \\\\\n&\\qquad at the split point,")
print("but this would significantly increase the equation count and may")
print("fragment the mathematical structure.")
print("\nThe current state (with \\\\ newlines added) is likely the best")
print("balance between raw-source readability and mathematical coherence.")
