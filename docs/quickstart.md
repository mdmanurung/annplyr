# Quickstart: from AnnData to a cohort summary

This page follows a small but realistic single-cell task: retain QC-passing
cells, attach experimental metadata, derive a marker score, summarize the
cohort, and prepare a long table for plotting. The examples execute during the
documentation build against a four-cell PBMC-like AnnData fixture.

## 1. Filter cells with metadata predicates

Importing `annplyr` registers the accessor. Multiple predicates in a tuple are
combined with `AND`, so the filtering criteria read like the QC rule itself.

```{testcode}
import annplyr as ap

qc = adata.ap.filter(
    obs=(
        ap.col("n_counts") >= 1_000,
        ap.col("pct_counts_mt") < 10,
    )
)

assert list(qc.obs_names) == ["cell_0", "cell_2", "cell_3"]
assert not qc.is_view
```

Axis-changing verbs return an independent AnnData by default. The original
object is still available for alternative QC thresholds or diagnostics.

## 2. Attach a sample sheet

AnnData-safe joins enrich metadata without manufacturing matrix rows. Declaring
the expected relationship turns a duplicated sample-sheet key into a clear
error instead of duplicated cells.

```{testcode}
sample_sheet = pd.DataFrame(
    {
        "batch": ["A", "B"],
        "condition": ["control", "stimulated"],
    }
)

annotated = qc.ap.left_join(
    sample_sheet,
    by="batch",
    relationship="many-to-one",
)

assert annotated.obs["condition"].tolist() == [
    "control",
    "control",
    "stimulated",
]
```

## 3. Derive metadata from expression values

`mutate()` writes metadata while matrix arguments remain read-only sources.
Here the two-gene score becomes an `obs` column; `X` itself is unchanged.

```{testcode}
scored = annotated.ap.mutate(
    x={"B_cell_score": (ap.col("MS4A1") + ap.col("CD79A")) / 2},
)

assert scored.obs["B_cell_score"].tolist() == [1.0, 3.0, 2.0]
assert "B_cell_score" not in annotated.obs
```

The same pattern works with a named layer by passing `layer="counts"`, or with
`raw`, `obsm`, and `varm` as explicit expression sources.

## 4. Summarize biological and experimental groups

One-off summaries accept `by=` directly. Metadata and matrix reductions can be
evaluated in the same call.

```{testcode}
cohort_summary = scored.ap.summarize(
    obs={
        "cells": ap.n(),
        "mean_counts": ap.mean("n_counts"),
    },
    x={"mean_MS4A1": ap.mean("MS4A1")},
    by=["condition", "cell_type"],
)

assert cohort_summary["cells"].tolist() == [2, 1]
assert cohort_summary["condition"].tolist() == ["control", "stimulated"]
```

Use persistent grouping when several steps share the same groups:

```{testcode}
ranked = (
    scored.ap.group_by(obs="condition")
    .mutate(obs={"within_condition": ap.min_rank("n_counts", descending=True)})
    .ungroup()
)

assert ranked.obs["within_condition"].tolist() == [2, 1, 1]
```

## 5. Prepare a plot-ready expression table

Select only the features required by the figure. `max_matrix_values=` bounds
the cumulative matrix projection before any source is read.

```{testcode}
plot_data = scored.ap.to_tidy(
    obs=["condition", "cell_type"],
    x=["MS4A1", "CD79A"],
    max_matrix_values=2 * scored.n_obs,
)

assert plot_data.shape == (6, 5)
assert list(plot_data.columns) == [
    "obs_name",
    "feature",
    "value",
    "condition",
    "cell_type",
]
```

`plot_data` is an ordinary pandas DataFrame ready for seaborn, plotnine,
Altair, statistical modelling, or export. `scored` remains a fully aligned
AnnData object for the rest of the analysis.

Next, follow the {doc}`tutorials`, learn the storage model in
{doc}`user_guide/concepts`, or look up an exact signature in {doc}`api`.
