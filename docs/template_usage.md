# Template Usage

`annplyr` is scaffolded in the style of
[`scverse/cookiecutter-scverse`](https://github.com/scverse/cookiecutter-scverse).
The template provides the baseline expectations for a public scverse ecosystem
package.

## Online Services

Enable these services for full infrastructure support:

- GitHub Actions for tests, docs, package build, and release workflows.
- GitHub Pages for hosted documentation, with pull request docs builds in
  GitHub Actions.
- Codecov for project and patch coverage reporting.
- pre-commit.ci for automatic hook checks on pull requests.
- PyPI trusted publishing for release uploads from GitHub Actions.

## Template Sync

The repository includes `.cruft.json` so future template updates can be checked
with cruft. Template updates should be reviewed as normal pull requests because
`annplyr` intentionally customizes the API shape away from scanpy-style
`pp`/`tl`/`pl` modules.

## Local Commands

```bash
pytest -q
uvx hatch run docs:build
prek run --all-files
uv build
uvx twine check --strict dist/*
```
