# annplyr

`annplyr` provides dataframe-style wrangling for AnnData objects through an
`adata.ap` accessor.

```python
import annplyr as ap

filtered = adata.ap.filter(obs=ap.col("batch") == "A", x=ap.col("GeneA") > 0)
plot_data = filtered.ap.to_df(obs=["batch"], x=["GeneA"])
```

## Development

This repository follows the scverse cookiecutter style for local tooling:

```bash
pytest
python -m build
```

The package metadata is managed through `pyproject.toml`, with Hatch
environments for tests and documentation.
