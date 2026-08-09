# External typing consumer

This is a standalone downstream project: it depends only on the installed
`annplyr` distribution and imports no private modules. The type-check job runs
its strict mypy configuration after the internal package check.

The narrow `anndata.*` override acknowledges that current stable AnnData does
not publish a PEP 561 marker. It does not ignore missing or untyped imports for
annplyr, pandas, or Narwhals, so an absent `py.typed` marker or a broken public
annotation fails this fixture. `typing.assert_type()` freezes accessor,
grouped, expression, join, extraction, and pipe inference.
