"""Install the bundled annplyr Agent Skill for Claude Code and Codex."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Literal

SKILL_NAME = "annplyr"
AgentName = Literal["claude", "codex"]
DEFAULT_AGENTS: tuple[AgentName, ...] = ("claude", "codex")


def bundled_skill_dir() -> Path:
    """Return the path to the Agent Skill bundled inside the package."""
    return Path(__file__).resolve().parent / "data"


def default_dest(agent: AgentName) -> Path:
    """Return the default personal skill destination for an agent."""
    if agent == "claude":
        return Path.home() / ".claude" / "skills" / SKILL_NAME
    if agent == "codex":
        codex_home = Path(os.environ.get("CODEX_HOME", Path.home() / ".codex"))
        return codex_home / "skills" / SKILL_NAME
    raise ValueError(f"Unsupported agent: {agent}")


def install_skill(dest: Path | None = None, *, agent: AgentName = "claude", force: bool = False) -> Path:
    """Copy the bundled annplyr skill into a Claude Code or Codex skills directory."""
    src = bundled_skill_dir()
    if not (src / "SKILL.md").is_file():
        raise FileNotFoundError(
            f"Bundled skill not found at {src}. The package may be installed without its skill data."
        )

    if dest is None:
        dest = default_dest(agent)

    if dest.exists():
        if not force:
            raise FileExistsError(f"{dest} already exists. Re-run with --force to overwrite.")
        shutil.rmtree(dest)

    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest, ignore=shutil.ignore_patterns(".ipynb_checkpoints", "__pycache__"))
    return dest


def install_skills(
    agents: Iterable[AgentName] = DEFAULT_AGENTS,
    *,
    dest: Path | None = None,
    force: bool = False,
) -> dict[AgentName, Path]:
    """Install the bundled skill for one or more supported agents."""
    selected = tuple(dict.fromkeys(agents))
    if dest is not None and len(selected) != 1:
        raise ValueError("--dest can only be used when exactly one --agent is selected.")

    return {agent: install_skill(dest=dest, agent=agent, force=force) for agent in selected}


def _display_agent(agent: AgentName) -> str:
    return {"claude": "Claude Code", "codex": "Codex"}[agent]


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``annplyr-install-skills``."""
    parser = argparse.ArgumentParser(
        prog="annplyr-install-skills",
        description="Install the annplyr Agent Skill for Claude Code and Codex.",
    )
    parser.add_argument(
        "--agent",
        action="append",
        choices=DEFAULT_AGENTS,
        help="Agent to install for. May be passed more than once. Defaults to Claude Code and Codex.",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=None,
        help="Destination directory. Use only with one --agent.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite an existing installation.")
    parser.add_argument(
        "--print-path",
        action="store_true",
        help="Print the bundled skill directory and exit without installing.",
    )
    args = parser.parse_args(argv)

    if args.print_path:
        print(bundled_skill_dir())
        return 0

    agents = tuple(args.agent) if args.agent else DEFAULT_AGENTS

    try:
        installed = install_skills(agents, dest=args.dest, force=args.force)
    except (FileExistsError, FileNotFoundError, ValueError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    for agent, dest in installed.items():
        print(f"Installed annplyr skill for {_display_agent(agent)} to {dest}")
    print("Restart the agent session for the installed skill to be discovered.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
