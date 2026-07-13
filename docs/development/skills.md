# Agent Skills

`annplyr` ships a bundled Agent Skill for Claude Code and Codex. The skill
teaches agents the `adata.ap` accessor, AnnData alignment rules, safe source
usage, plot-ready extraction, and sparse/backed-data constraints.

Install both Claude Code and Codex copies:

```bash
annplyr-install-skills
```

Install only the Claude Code copy:

```bash
annplyr-install-skills --agent claude
```

Install only the Codex copy:

```bash
annplyr-install-skills --agent codex
```

Refresh an existing copy after upgrading `annplyr`:

```bash
annplyr-install-skills --force
```

Print the bundled skill path without installing it:

```bash
annplyr-install-skills --print-path
```

The skill uses a compact router file plus reference pages that are loaded on
demand, so it stays small until an agent needs specific API or safety details.
