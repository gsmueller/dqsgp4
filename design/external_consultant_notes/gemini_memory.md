# Gemini Tool Usage Memory & Advice

This document serves as a persistent cache of operational wisdom for tool interaction within this workspace.

## 🚀 General Efficiency Strategies
- **Minimize Turns:** Each turn is a permanent addition to session history. Combine search and read operations in a single turn.
- **Strategic Grepping:** Use `grep_search` with `context`, `before`, or `after` parameters to identify `old_string` targets immediately, bypassing the need for a follow-up `read_file`.
- **Parallelism:** Dispatch independent tool calls (e.g., searching multiple directories or reading multiple files) in parallel. Set `wait_for_previous: true` ONLY for sequential dependencies.
- **Surgical Reads:** Use `start_line` and `end_line` aggressively. Reading more than necessary wastes context and risks hitting the 2000-line truncation limit.

## 🛠️ Tool-Specific Advice

### `grep_search` (Primary Discovery)
- **Scope Limitation:** Use `include_pattern` (e.g., `src/**/*.cpp`) to reduce noise in large codebases.
- **Regex Precision:** Use `\bSymbolName\b` for exact matches to avoid partial string hits.
- **Efficiency:** Favor this over `run_shell_command("grep ...")` for faster performance and automatic output limiting.

### `read_file` (Deep Dive)
- **Incremental Reading:** If a file is large and `grep_search` failed, read the first 100 lines and look for a table of contents or structural markers.

### `replace` (Surgical Updates)
- **Unambiguous Context:** Provide 3–5 lines of surrounding code in `old_string`. If `replace` fails, it's usually because the context was too narrow or shared by multiple blocks.
- **One Edit Per File:** NEVER attempt multiple `replace` calls on the same file in one turn; this causes race conditions.
- **Idiomatic Completeness:** Ensure the `new_string` follows local formatting (e.g., indentation, braces) exactly.

### `run_shell_command` (System Interaction)
- **PowerShell Limitations:** `cat << 'EOF'` (heredocs) fail in the Windows PowerShell environment. Use `write_file` to create multi-line scripts or configuration files instead.
- **Sequential Execution:** Since `&&` is not supported, use `;` to separate commands or dispatch them in separate tool calls with `wait_for_previous: true`.
- **Quiet Mode:** Always prefer silent flags (e.g., `npm install --silent`, `git --no-pager`) to minimize context bloat.
- **Non-Interactive:** Ensure commands do not prompt for input. Use ` -y` or `force` flags where available.
- **Explanation Rule:** Before running commands that modify the state (filesystem, git, etc.), provide a concise explanation of the intent.

## 🧪 Numerical Verification Pattern
- **Python/NumPy:** For complex mathematical verification (FFT, Kepler solving, integration), write a temporary Python script using `write_file` and execute it with `run_shell_command`.
- **Environment:** `numpy` is available and should be used for vectorized orbital calculations.
- **Cleanup:** Always `rm` temporary verification scripts after use to keep the workspace clean.

### `invoke_agent` (Delegation)
- **Context Compression:** Delegate repetitive batch tasks (e.g., "Fix lint errors in 10 files") or speculative research. The subagent's execution history is collapsed, keeping the main loop lean.

### `save_memory` (Persistence)
- **Project Facts:** Use `scope: "project"` for facts about the codebase architecture (e.g., "The $\eta$-denominator policy is $\eta^{2p-3}$").
- **Global Preferences:** Use `scope: "global"` for general peer-programming preferences.

## ⚖️ Engineering Mandates
- **Validation is Finality:** A task is incomplete until the behavioral correctness of the change is verified (tests) and its structural integrity is confirmed (lint/type-check).
- **Research First:** Never act on a bug report without first empirically reproducing the failure.
- **Contextual Precedence:** Instructions in `GEMINI.md` take absolute precedence over all defaults.
