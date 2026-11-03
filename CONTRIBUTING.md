# Contributing Guide

This project follows the scverse cookiecutter development style: Hatch
environments, uv-based installs, pytest, Ruff, prek/pre-commit, GitHub Actions,
GitHub Pages, and Codecov.

## Development Environment

Use Hatch through `uvx` when you want the same environments as CI:

```bash
uvx hatch env create hatch-test
uvx hatch test
uvx hatch run docs:build
```

Use uv when you want a single local environment:

```bash
uv sync --all-groups
uv run pytest
```

## Code Style

Run the repository hooks before committing:

```bash
prek run --all-files
```

The hooks are defined in `.pre-commit-config.yaml`; Ruff configuration lives in
`pyproject.toml`.

## Tests

Every behavior change should include tests. For annplyr, this usually means:

- metadata tests for `obs` and `var` behavior,
- matrix tests for dense, sparse, and layer-backed sources,
- alignment tests that verify AnnData containers keep their expected shapes,
- error tests for invalid selectors, sources, axes, and ambiguous names.

Fast local gate:

```bash
pytest -q
```

Release-quality gate:

```bash
uvx hatch test --all
```

## Documentation

Documentation is built with Sphinx and MyST-NB:

```bash
uvx hatch run docs:build
```

The seven pages in `docs/*.md` are MyST-NB notebooks. The build executes every
cell, so each example and figure is checked, and a failing cell fails the build.
Hidden `remove-cell` cells assert that outputs match. Shared setup is in
`docs/_setup.py`. The build downloads PBMC3k into `data/` on first use.

On pushes to `main`, the `Docs` workflow uploads `docs/_build/html` and deploys
the Sphinx site to GitHub Pages.

## Releasing

Releases are published from GitHub releases through the `release.yaml` workflow
and PyPI trusted publishing.

Before creating a release:

1. Update the version in `pyproject.toml`.
2. Update `CHANGELOG.md` and `CITATION.cff` with the actual publication date.
3. Run tests, docs, hooks, and package checks.
4. Create a GitHub release tagged as `vX.Y.Z`.

PyPI must have a trusted publisher configured for repository
`mdmanurung/annplyr`, workflow `release.yaml`, and environment `pypi`.

After publishing, confirm that the release workflow succeeds, install the exact
version from PyPI in a clean environment, run a minimal AnnData smoke workflow,
and verify the published documentation and release links.
