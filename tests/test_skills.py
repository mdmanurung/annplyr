from __future__ import annotations

import tomllib
from pathlib import Path

import pytest

from annplyr._skills import install


def test_bundled_skill_directory_contains_router_and_references() -> None:
    skill_dir = install.bundled_skill_dir()
    skill_text = (skill_dir / "SKILL.md").read_text()

    assert (skill_dir / "SKILL.md").is_file()
    assert (skill_dir / "agents" / "openai.yaml").is_file()
    assert (skill_dir / "references" / "quickstart.md").is_file()
    assert "name: annplyr" in skill_text
    assert "license: BSD-3-Clause" in skill_text


def test_skill_trigger_description_covers_core_user_phrases() -> None:
    skill_text = (install.bundled_skill_dir() / "SKILL.md").read_text()
    description = next(line for line in skill_text.splitlines() if line.startswith("description: "))

    for phrase in ["adata.ap", "to_df", "to_tidy", "plotnine", "single-cell"]:
        assert phrase in description


def test_skill_examples_only_use_supported_annplyr_expression_methods() -> None:
    quickstart = (install.bundled_skill_dir() / "references" / "quickstart.md").read_text()

    assert ".log1p()" not in quickstart


def test_install_skill_copies_bundle_and_refuses_existing_destination(tmp_path: Path) -> None:
    dest = tmp_path / "skills" / "annplyr"

    installed = install.install_skill(dest=dest)

    assert installed == dest
    assert (dest / "SKILL.md").is_file()
    assert (dest / "references" / "api-patterns.md").is_file()
    with pytest.raises(FileExistsError, match="--force"):
        install.install_skill(dest=dest)

    sentinel = dest / "old.txt"
    sentinel.write_text("stale")
    install.install_skill(dest=dest, force=True)
    assert not sentinel.exists()


def test_agent_destinations_use_expected_personal_skill_roots(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("CODEX_HOME", str(tmp_path / "codex-home"))

    assert install.default_dest("claude") == tmp_path / "home" / ".claude" / "skills" / "annplyr"
    assert install.default_dest("codex") == tmp_path / "codex-home" / "skills" / "annplyr"


def test_main_installs_selected_agent_and_prints_bundled_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("CODEX_HOME", raising=False)

    assert install.main(["--print-path"]) == 0
    assert capsys.readouterr().out.strip() == str(install.bundled_skill_dir())

    assert install.main(["--agent", "codex"]) == 0
    output = capsys.readouterr().out

    assert "Installed annplyr skill for Codex" in output
    assert (tmp_path / "home" / ".codex" / "skills" / "annplyr" / "SKILL.md").is_file()
    assert not (tmp_path / "home" / ".claude" / "skills" / "annplyr").exists()


def test_console_script_is_declared() -> None:
    pyproject = tomllib.loads(Path("pyproject.toml").read_text())

    assert pyproject["project"]["scripts"]["annplyr-install-skills"] == "annplyr._skills.install:main"
