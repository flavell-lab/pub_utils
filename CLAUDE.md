# CLAUDE.md - CRITICAL RULES
This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
IMPORTANT: You MUST adhere to the Token Efficiency Guidelines below. 
Prioritize these constraints even if they conflict with your default verbosity.


## Claude Code Session Guidelines

1. **Daily logging**: Keep track of output files and the logic behind their creation. Save logs under `claudecode/yyyy-mm-dd.py` (e.g., `claudecode/2026-01-14.py`). Prepare the daily log file when I have <3% context left until auto-compact.

2. **Ask, don't guess**: When requirements are ambiguous or in conflict with each other, ask the user for clarification rather than interpolating or guessing.

3. **Token efficiency**: Use simple queries for easy tasks to save tokens. Avoid over-exploring when the task is straightforward.

4. **Clear documentation**: The logic behind the data ETL and connectome assembly needs to be extremely clearly documented for future contributors to reproduce the outputs.

5. **Continuity between sessions**: At the beginning of each claude session, load in the last log file from `claudecode` so that you know where to pick up from.

6. **Error handling**: Minimize try/except patterns in `src/` functions. If a try block is necessary, emit a `warnings.warn()` message on failure rather than silently continuing. This ensures users are aware when something unexpected happens.


## Token Efficiency Guidelines

### 1. Context Control
- **Minimize file additions:** Only add files strictly necessary for the current task using `+`. 
- **Remove files when done:** Use `-` to remove files from context once they are no longer being modified to keep the context window lean.
- **Scope by Subdirectory:** If working on a massive repo, start the `claude` CLI within a specific subdirectory to limit the initial file tree indexing.

### 2. Prompting & Output Discipline
- **Task Specificity:** Avoid broad requests like "refactor this." Use "refactor the `getUsers` method in `api.ts` to use the new wrapper."
- **Sequential Workflow:** Complete one feature and verify it before moving to the next. Long "todo lists" in a single session explode token usage.
- **Request Conciseness:** Use "be concise" or "skip explanations" if you only need the code. This significantly reduces completion tokens.
- **Prefer Diffs:** Encourage Claude to show snippets or diffs for large files instead of rewriting the entire file in the chat.

### 3. Session Management
- **Use `/clear` for pivot points:** When switching from a frontend task to a backend task, use `/clear`. This wipes the history and resets the token count for a fresh, high-performance start.
- **Strategic Compaction:** Only use `/compact` when necessary, and always specify what to preserve (e.g., `/compact keep the database schema details`).
- **Terminal Awareness:** Avoid running commands that produce massive stdout (like verbose test logs). Filter or pipe the output so Claude doesn't "read" thousands of lines of logs.

### 4. Memory Management
- **Leverage CLAUDE.md:** Keep core architectural rules and frequent patterns here. It’s more efficient for Claude to read this once than to have you repeat instructions in every session.