# A sample-aware marker workflow

Single-cell wrangling often crosses several kinds of data at once. Cell-level
QC lives in `obs`, marker values live in `X` or a layer, experimental condition
lives in a separate sample sheet, and the final figure needs a plain table.
This vignette keeps those steps in one readable workflow without detaching the
metadata from its AnnData object prematurely.

The documentation fixture represents four PBMC-like cells, three marker genes,
two donors, and a raw-count layer. The values are synthetic so the example is
small and deterministic; the structure mirrors an ordinary annotated cohort.

## Start by checking the inputs

Before selecting marker genes, distinguish missing features from genes that are
present under a different spelling. `feature_present()` reports the complete
request without changing the object.

```{testcode}
import annplyr as ap

markers = ["MS4A1", "CD79A", "NKG7", "IL7R"]
presence = ap.feature_present(adata, markers)

assert presence.found_features == ["MS4A1", "CD79A", "NKG7"]
assert presence.missing_features == ["IL7R"]
```

This is a useful validation boundary for a configurable pipeline: report the
missing genes once, then pass only the confirmed features to downstream verbs.

## Attach donor-level metadata safely

With the feature set validated, experimental annotations usually arrive as one
row per sample or donor. A left join adds those columns to `obs` while keeping
every cell and every aligned matrix row in place.

```{testcode}
sample_sheet = pd.DataFrame(
    {
        "donor_id": ["donor_1", "donor_2"],
        "condition": ["control", "stimulated"],
        "collection_day": [0, 7],
    }
)

cohort = adata.ap.left_join(
    sample_sheet,
    by="donor_id",
    relationship="many-to-one",
)

assert cohort.n_obs == adata.n_obs
assert cohort.obs["condition"].tolist() == [
    "control",
    "stimulated",
    "control",
    "stimulated",
]
assert "condition" not in adata.obs
```

Declaring the join relationship is more than documentation. If the sample
sheet accidentally contains two rows for one donor, annplyr raises
`JoinRelationshipError` instead of duplicating cells that have no corresponding
rows in `X`, layers, or embeddings.

## Apply the cell-level QC rule

QC thresholds remain ordinary expressions over `obs`. A tuple means that every
predicate must pass.

```{testcode}
qc = cohort.ap.filter(
    obs=(
        ap.col("n_counts") >= 1_000,
        ap.col("pct_counts_mt") < 10,
    )
)

assert list(qc.obs_names) == ["cell_0", "cell_2", "cell_3"]
assert qc.shape == (3, 3)
assert not qc.is_view
```

`filter()` returns an independent AnnData by default and subsets all aligned
containers through the same integer positions. That makes the ownership
boundary clear: `cohort` still represents the unfiltered data, while `qc` can
be passed to the next stage.

## Derive a marker score without rewriting the matrix

Matrix arguments to `mutate()` are read-only expression sources. On that
filtered object, the following call reads two selected genes and writes one new
`obs` column; it does not replace or densify `X`.

```{testcode}
scored = qc.ap.mutate(
    x={"B_cell_score": (ap.col("MS4A1") + ap.col("CD79A")) / 2},
    max_matrix_values=2 * qc.n_obs,
)

assert scored.obs["B_cell_score"].tolist() == [1.0, 3.0, 2.0]
assert np.array_equal(scored.X, qc.X)
```

For an AnnData object whose `X` stores normalized values and whose
`layers["counts"]` stores counts, add `layer="counts"` to evaluate the same
expression against that layer. `raw=`, keyed `obsm=`, and keyed `varm=` sources
make other aligned measurements available without inventing separate helper
APIs.

## Keep grouping explicit across related operations

With the score in place, persistent grouping is useful when several operations
share the same unit of comparison. Here it stays attached while cells are ranked
within condition. The workflow then removes that grouping explicitly before
computing one summary for each condition and cell type.

```{testcode}
grouped = scored.ap.group_by(obs="condition")
ranked = grouped.mutate(
    obs={"within_condition_rank": ap.min_rank("n_counts", descending=True)}
)

assert ranked.group_vars() == ["condition"]
assert ranked.ungroup().obs["within_condition_rank"].tolist() == [2, 1, 1]

cell_type_summary = ranked.ungroup().ap.summarize(
    obs={
        "cells": ap.n(),
        "mean_counts": ap.mean("n_counts"),
        "mean_B_cell_score": ap.mean("B_cell_score"),
    },
    x={"mean_MS4A1": ap.mean("MS4A1")},
    by=["condition", "cell_type"],
)

assert cell_type_summary["cells"].tolist() == [2, 1]
```

Grouped AnnData verbs keep returning a grouped wrapper, while summaries return
pandas DataFrames. Calling `ungroup()` marks the boundary where ordinary
AnnData chaining resumes.

## Build the analysis object and the plotting table

After ungrouping, the AnnData result can be narrowed to the metadata and genes
required by the next analysis. A separate long table then provides a natural
handoff to a plotting library.

```{testcode}
analysis = ranked.ungroup().ap.select(
    obs=[
        "donor_id",
        "condition",
        "cell_type",
        "n_counts",
        "B_cell_score",
        "within_condition_rank",
    ],
    x=["MS4A1", "CD79A", "NKG7"],
)

plot_data = analysis.ap.to_tidy(
    obs=["donor_id", "condition", "cell_type"],
    x=["MS4A1", "CD79A", "NKG7"],
    max_matrix_values=3 * analysis.n_obs,
)

assert analysis.shape == (3, 3)
assert plot_data.shape == (9, 6)
```

For example, seaborn can consume `plot_data` directly:

```python
import seaborn as sns

sns.catplot(
    data=plot_data,
    x="cell_type",
    y="value",
    hue="condition",
    col="feature",
    kind="strip",
    sharey=False,
)
```

The plotting dependency stays outside annplyr, and the table contains only the
requested nine matrix values. The finite budget documents that expectation and
causes an oversized request to fail before the first matrix read.

## Keep domain-specific steps composable

More generally, `pipe()` is the escape hatch for project-specific functions. A
function can accept AnnData, use ordinary Python or another scverse package, and
return any type appropriate for the next step.

```{testcode}
def donor_counts(data):
    return data.ap.count("donor_id", sort=True)


counts = analysis.ap.pipe(donor_counts)
assert counts["n"].sum() == analysis.n_obs
```

At this point `analysis` remains a standard AnnData object suitable for Scanpy,
CellRank, scvi-tools, serialization, or another `.ap` pipeline. annplyr owns the
wrangling semantics; the surrounding scientific workflow keeps ownership of
normalization, modelling, visualization, and biological interpretation.

## Complete pipeline

The central steps can be read as one compact workflow:

```python
analysis = (
    adata.ap.left_join(sample_sheet, by="donor_id", relationship="many-to-one")
    .ap.filter(
        obs=(
            ap.col("n_counts") >= 1_000,
            ap.col("pct_counts_mt") < 10,
        )
    )
    .ap.mutate(
        x={"B_cell_score": (ap.col("MS4A1") + ap.col("CD79A")) / 2},
        max_matrix_values=2 * adata.n_obs,
    )
    .ap.select(
        obs=["donor_id", "condition", "cell_type", "B_cell_score"],
        x=["MS4A1", "CD79A", "NKG7"],
    )
)
```

The individual calls remain ordinary, inspectable operations, but the complete
pipeline communicates the analysis more directly than a sequence of manual
boolean masks, `.obs` assignments, pandas merges, and synchronized AnnData
slices.
