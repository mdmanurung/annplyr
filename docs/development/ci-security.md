# CI security policy

The GitHub Actions surface is intentionally small and immutable. Every external
action is pinned to a full 40-character commit SHA, with the corresponding
upstream release retained as an inline comment. The release tags below were
resolved through the GitHub Git references API on 2026-08-09; annotated tags
were peeled to their commit objects before acceptance.

## Accepted action pins

| Action | Upstream release | Immutable commit | annplyr use |
|---|---|---|---|
| [`actions/checkout`](https://github.com/actions/checkout/releases/tag/v7.0.1) | `v7.0.1` | `3d3c42e5aac5ba805825da76410c181273ba90b1` | all jobs that read the repository |
| [`astral-sh/setup-uv`](https://github.com/astral-sh/setup-uv/releases/tag/v9.0.0) | `v9.0.0` | `c771a70e6277c0a99b617c7a806ffedaca235ff9` | Python, uv, Hatch, docs, tests, and builds |
| [`actions/upload-artifact`](https://github.com/actions/upload-artifact/releases/tag/v7.0.1) | `v7.0.1` | `043fb46d1a93c77aae656e7c1c64a875d1fc6a0a` | integration failure diagnostics |
| [`actions/configure-pages`](https://github.com/actions/configure-pages/releases/tag/v6.0.0) | `v6.0.0` | `45bfe0192ca1faeb007ade9deae92b16b8254a0d` | Pages configuration |
| [`actions/upload-pages-artifact`](https://github.com/actions/upload-pages-artifact/releases/tag/v5.0.0) | `v5.0.0` | `fc324d3547104276b827a68afc52ff2a11cc49c9` | rendered Pages artifact |
| [`actions/deploy-pages`](https://github.com/actions/deploy-pages/releases/tag/v5.0.0) | `v5.0.0` | `cd2ce8fcbc39b97be8ca5fce6e763baed58fa128` | main-branch Pages deployment |
| [`pypa/gh-action-pypi-publish`](https://github.com/pypa/gh-action-pypi-publish/releases/tag/v1.14.2) | `v1.14.2` | `dc37677b2e1c63e2034f94d8a5b11f265b73ba33` | PyPI trusted publishing |
| [`codecov/codecov-action`](https://github.com/codecov/codecov-action/releases/tag/v7.0.0) | `v7.0.0` | `fb8b3582c8e4def4969c97caa2f19720cb33a72f` | advisory OIDC coverage upload |
| [`re-actors/alls-green`](https://github.com/re-actors/alls-green/releases/tag/v1.2.2) | `v1.2.2` | `05ac9388f0aebcb5727afa17fcccfecd6f8ec5fe` | aggregate stable and advisory test jobs |
| [`zizmorcore/zizmor-action`](https://github.com/zizmorcore/zizmor-action/releases/tag/v0.6.2) | `v0.6.2` | `3dc1ecc9bcb9e94e9b2c709687979e1298497054` | blocking workflow audit with zizmor 1.29.0 |

The accepted major updates retain the inputs annplyr uses. `setup-uv` v9 changes
only the default cache-pruning policy relevant here; release jobs disable that
cache entirely. The Pages actions target GitHub-hosted runners and GitHub.com,
and the artifact jobs continue to use their default archived upload path.

## Permissions and credentials

Build, typing, integration, documentation, benchmark, and test-discovery jobs
receive only `contents: read`. The test matrix adds `id-token: write` solely for
Codecov OIDC. The Pages deployment adds only `pages: write` and
`id-token: write`; the release job adds only `id-token: write` for PyPI trusted
publishing. The release workflow otherwise defaults to no token permissions.

Every checkout sets `persist-credentials: false`. Release builds disable the uv
cache, and release concurrency prevents simultaneous trusted-publishing jobs.
Matrix-derived Hatch environment names cross into shell steps through a quoted
environment variable rather than direct template expansion.

## Continuous audit and dependency updates

`.github/workflows/security.yaml` runs the pedantic zizmor persona on pull
requests, pushes to `main`, a weekly schedule, and manual dispatches. It pins
both `zizmor-action` and zizmor itself, fails on findings, and uses console plus
bounded annotation output instead of requesting `security-events: write`.
Online audits use the job's read-only GitHub token.

The pre-change offline audit with zizmor 1.29.0 reported 60 findings: 27 high,
9 medium, 19 low, and 5 informational. Pinning every action, disabling checkout
credential persistence, constraining permissions, documenting write scopes,
removing direct matrix interpolation, disabling the release cache, and naming
jobs reduced the same pedantic audit to zero findings.

Dependabot monitors only the `github-actions` ecosystem. It opens one grouped
weekly update after a seven-day cooldown, with at most five open update pull
requests. It does not monitor Python or pre-commit dependencies, so it does not
overlap the repository's existing dependency-update mechanisms.

## Deliberately excluded template changes

Closed PR #8 was reviewed only as source material. This focused replacement
does not adopt its package-version and dependency resets, one-percent Codecov
target, documentation/template churn, pre-commit auto-fix hook, or removal of
strict Twine validation. The 0.3.0 package metadata, 85% branch-coverage floor,
release-on-published trigger, Pages deployment, and PyPI OIDC flow remain
unchanged.
