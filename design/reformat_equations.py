"""
Reformat equations in the Brouwer-Hori 1961 markdown.
Adds line endings after each \\\\ continuation marker within $$ equation blocks.
Does not add repeats (checks if newline already follows \\\\).
"""

import re
import sys

input_path = r"sgp4_references\vallado_celestrak\documentation\SGP4\Brouwer_Hori 1961\Brouwer_Hori_1961a_Atmospheric_Drag_AJ66_193.md"

with open(input_path, 'r', encoding='utf-8') as f:
    content = f.read()

# Strategy: find all $$ ... $$ blocks and within them,
# replace \\ followed by non-newline with \\\n
# But we need to be careful not to break \\[ or \\\\ or other LaTeX sequences.

# The pattern: \\ followed by a space or &, which indicates a line continuation
# in an aligned environment. We want to add a newline after the \\ if one isn't there.

# More precisely: in LaTeX aligned environments, \\ marks a line break.
# We want the raw markdown to have an actual newline after each \\.
# Pattern: \\ followed by optional whitespace then & or \end or ]
# But simplest: just add newline after \\ when it's not already followed by newline.

# Find $$ blocks
def reformat_equation(match):
    """Add newlines after \\\\ in equation blocks."""
    text = match.group(0)

    # Only process blocks that contain \begin{aligned}
    # (single-line equations don't need this)
    if r'\begin{aligned}' not in text:
        return text

    # Replace \\ followed by non-newline with \\\n
    # But preserve \\ at end of line (already has newline)
    # Pattern: \\ followed by space(s) then non-newline content
    result = re.sub(r'\\\\(\s*)(?!\n)', r'\\\\\n', text)

    return result

# Match $$ ... $$ blocks (including multiline)
# Use DOTALL so . matches newlines
pattern = re.compile(r'\$\$.*?\$\$', re.DOTALL)

new_content = pattern.sub(reformat_equation, content)

# Count changes
original_lines = content.count('\n')
new_lines = new_content.count('\n')
print(f"Original: {original_lines} lines")
print(f"New:      {new_lines} lines")
print(f"Added:    {new_lines - original_lines} line breaks")

# Write output
with open(input_path, 'w', encoding='utf-8') as f:
    f.write(new_content)

print(f"\nFile updated: {input_path}")
