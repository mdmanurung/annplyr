# annplyr

`annplyr` provides tidy, dataframe-style wrangling for AnnData objects through
an `adata.ap` accessor. It is built for single-cell workflows where metadata,
expression matrices, layers, embeddings, and loadings need to be queried
together without losing AnnData alignment.

```python
import annplyr as ap

filtered = adata.ap.filter(
    obs=ap.col("batch") == "A",
    x=ap.col("MS4A1") > 0,
)

plot_data = filtered.ap.to_tidy(obs=["cell_type"], x=["MS4A1", "CD79A"])
```

## Start Here

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} Install annplyr
:link: installation
:link-type: doc

Set up the package and development environment.
:::

:::{grid-item-card} Quickstart
:link: quickstart
:link-type: doc

Run the first AnnData wrangling examples.
:::

:::{grid-item-card} Tutorials
:link: tutorials
:link-type: doc

Work through the notebooks, including the comprehensive tour.
:::

:::{grid-item-card} API Reference
:link: api
:link-type: doc

Look up accessor methods, expressions, selectors, helpers, and errors.
:::

::::

```{toctree}
:caption: Get Started
:maxdepth: 1

installation
quickstart
tutorials
```

```{toctree}
:caption: User Guide
:maxdepth: 1

user_guide/concepts
user_guide/verbs
user_guide/grouping
user_guide/joins
user_guide/tidy_tables
user_guide/scverse_safety
user_guide/utilities
```

```{toctree}
:caption: API Reference
:maxdepth: 1

api
```

```{toctree}
:caption: Development
:maxdepth: 1

development/skills
contributing
template_usage
```

```{toctree}
:caption: Project
:maxdepth: 1

roadmap
changelog
references
```
