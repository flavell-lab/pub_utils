---
description: Logs progress to a daily file and prepares for a /clear reset.
---

# Action: Document & Checkpoint

1. **Log to File:** Use your tools to ensure a directory named `./claudecode/` exists. 
2. **Write Summary:** Create or append to a file named `./claudecode/{{DATE}}.md` (use today's date in YYYY-MM-DD format). 
   - Write a concise log of what was changed, which files were touched, and the current build status.
   - If the file exists, append a new "Session Checkpoint" entry with a timestamp.

3. **Provide CLI Response:** Once the file is written, output the following to the terminal:

### 📝 Task Documented
**Log Location:** `./claudecode/{{DATE}}.md`

**Next Step Prompt:**
> "I am resuming work from a previous session. Last progress: [Paste summary of current state]. Please re-index these files: [List of files to @mention]."

**Context to Restore:** (List the files I should add back with `+` or `@` when I restart)

---
**Summary complete and logged.** You can now run `/clear` to reset your token usage. Copy the 'Next Step Prompt' above to resume.