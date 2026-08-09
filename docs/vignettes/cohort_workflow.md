# A sample-aware marker workflow

Single-cell wrangling crosses several kinds of data at once. Cell-level QC lives
in `obs`, marker values live in `X` or `raw`, experimental condition lives in a
separate sample sheet, and the final figure needs a plain table. This vignette
runs all of that on Scanpy's PBMC3K without detaching metadata from its AnnData
object until the very last step.

```{note}
PBMC3K is a **single donor**. To demonstrate sample-aware wrangling, this page
assigns three synthetic donor labels by cell position and treats two of them as
stimulated. Those two columns — `donor_id` and `condition` — are the only
invented values on the page; every count, expression value, and cell-type label
below is real.
```

## Start by validating the panel

Before selecting marker genes, separate genuinely missing features from ones
present under a different spelling. `feature_present()` reports the whole
request without touching the object:

```{testcode}
markers = ["MS4A1", "CD79A", "NKG7", "LYZ", "CD4"]

print(ap.feature_present(adata, markers).to_frame())
```

```{testoutput}
  feature  found wrong_case_found
0   MS4A1   True             <NA>
1   CD79A   True             <NA>
2    NKG7   True             <NA>
3     LYZ  False             <NA>
4     CD4  False             <NA>
```

Two of five markers are missing — because `feature_present()` checks
`var_names`, and PBMC3K's `X` holds only the 1,838 highly variable genes.
Validate against the source you actually intend to read:

```{testcode}
print(ap.feature_present(adata.raw.to_adata(), markers).to_frame())
```

```{testoutput}
  feature  found wrong_case_found
0   MS4A1   True             <NA>
1   CD79A   True             <NA>
2    NKG7   True             <NA>
3     LYZ   True             <NA>
4     CD4   True             <NA>
```

All five are in `.raw`. That is the validation boundary for a configurable
pipeline: report once, then pass only confirmed features downstream.

## Attach donor-level metadata safely

Experimental annotations arrive as one row per donor. A left join adds those
columns to `obs` while keeping every cell and every aligned matrix row in place:

```{testcode}
import pandas as pd

adata.obs["donor_id"] = [f"donor_{i % 3 + 1}" for i in range(adata.n_obs)]

sample_sheet = pd.DataFrame(
    {
        "donor_id": ["donor_1", "donor_2", "donor_3"],
        "condition": ["control", "stimulated", "stimulated"],
        "collection_day": [0, 7, 7],
    }
)

cohort = adata.ap.left_join(sample_sheet, by="donor_id", relationship="many-to-one")

print(cohort.obs[["louvain", "donor_id", "condition", "collection_day"]].head(3))
```

```{testoutput}
                      louvain donor_id   condition  collection_day
index                                                             
AAACATACAACCAC-1  CD4 T cells  donor_1     control               0
AAACATTGAGCTAC-1      B cells  donor_2  stimulated               7
AAACATTGATCAGC-1  CD4 T cells  donor_3  stimulated               7
```

Declaring `relationship="many-to-one"` is more than documentation. If the sample
sheet accidentally held two rows for one donor, annplyr raises
`JoinRelationshipError` rather than duplicating cells that have no corresponding
rows in `X`, layers, or embeddings.

## Apply the cell-level QC rule

QC thresholds stay ordinary expressions over `obs`; a tuple means every
predicate must pass:

```{testcode}
qc = cohort.ap.filter(
    obs=(
        ap.col("n_genes") >= 500,
        ap.col("percent_mito") < 0.03,
    )
)

print(qc.shape, qc.is_view)
```

```{testoutput}
(2117, 1838) False
```

`filter()` returns an independent AnnData and subsets every aligned container
through the same integer positions, so `cohort` still represents the unfiltered
data while `qc` moves to the next stage.

## Derive a marker score without rewriting the matrix

Matrix arguments to `mutate()` are read-only sources. This call reads two genes
out of the 13,714 in `.raw` and writes one `obs` column:

```{testcode}
scored = qc.ap.mutate(
    raw={"b_score": (ap.col("MS4A1") + ap.col("CD79A")) / 2},
    max_matrix_values=2 * qc.n_obs,
)

print(scored.obs[["donor_id", "condition", "louvain", "b_score"]].head(3).round(3))
```

```{testoutput}
                 donor_id   condition          louvain  b_score
index                                                          
AAACATTGATCAGC-1  donor_3  stimulated      CD4 T cells      0.0
AAACCGTGCTTCCG-1  donor_1     control  CD14+ Monocytes      0.0
AAACCGTGTATGCG-1  donor_2  stimulated         NK cells      0.0
```

Add `layer="counts"` to evaluate the same expression against a named layer.
`raw=`, keyed `obsm=`, and keyed `varm=` make other aligned measurements
available without a separate helper API for each.

## Keep grouping explicit across related operations

Persistent grouping pays off when several operations share a unit of
comparison. Here the ranking and the selection both happen within donor:

```{testcode}
ranked = scored.ap.group_by(obs="donor_id").mutate(
    obs={"depth_rank": ap.min_rank("n_counts", descending=True)}
)

deepest = ranked.filter(obs=ap.col("depth_rank") <= 2).ungroup()

print(deepest.obs[["donor_id", "n_counts", "depth_rank"]])
```

```{testoutput}
                 donor_id  n_counts  depth_rank
index                                          
CAGGTTGAGGATCT-1  donor_3    8011.0         1.0
CGATCAGATGTGAC-1  donor_3    6908.0         2.0
ACGAACTGGCTATG-1  donor_1    8875.0         1.0
CATACTTGGGTTAC-1  donor_1    7167.0         2.0
ACGAGGGACAGGAG-1  donor_2    7928.0         2.0
GGGCCAACCTTGGA-1  donor_2    8415.0         1.0
```

Two cells per donor, ranked within donor, with the original cell order
preserved. Grouped AnnData verbs keep returning a grouped wrapper; `ungroup()`
marks where ordinary AnnData chaining resumes.

## Summarize the cohort

One call crosses condition and cell type, and mixes metadata with expression:

```{testcode}
summary = ranked.ungroup().ap.summarize(
    obs={"cells": ap.n(), "mean_b_score": ap.mean("b_score")},
    raw={"mean_LYZ": ap.mean("LYZ")},
    by=["condition", "louvain"],
)

print(summary.round(3).head(8))
```

```{testoutput}
    condition            louvain  cells  mean_b_score  mean_LYZ
0  stimulated        CD4 T cells    648         0.030     0.442
1     control    CD14+ Monocytes    118         0.029     3.733
2  stimulated           NK cells     89         0.019     0.373
3  stimulated        CD8 T cells    159         0.027     0.387
4     control  FCGR3A+ Monocytes     44         0.039     2.108
5  stimulated            B cells    183         1.239     0.425
6     control        CD4 T cells    337         0.043     0.462
7  stimulated    CD14+ Monocytes    231         0.038     3.755
```

The B-cell score and monocyte `LYZ` track cell type, not condition — exactly
what should happen when the condition labels are arbitrary. A real contrast
would show up as a difference between the `control` and `stimulated` rows of the
same cell type.

## Build the analysis object and the plotting table

```{testcode}
analysis = ranked.ungroup().ap.select(
    obs=["donor_id", "condition", "louvain", "n_counts", "b_score", "depth_rank"],
    x=["MS4A1", "CD79A", "NKG7"],
)

plot_data = analysis.ap.to_tidy(
    obs=["donor_id", "condition", "louvain"],
    raw=["MS4A1", "CD79A", "NKG7"],
    max_matrix_values=3 * analysis.n_obs,
)

print(analysis.shape)
print(plot_data.head())
print(plot_data.shape)
```

```{testoutput}
(2117, 3)
           obs_name feature  value donor_id   condition          louvain
0  AAACATTGATCAGC-1   MS4A1    0.0  donor_3  stimulated      CD4 T cells
1  AAACATTGATCAGC-1   CD79A    0.0  donor_3  stimulated      CD4 T cells
2  AAACATTGATCAGC-1    NKG7    0.0  donor_3  stimulated      CD4 T cells
3  AAACCGTGCTTCCG-1   MS4A1    0.0  donor_1     control  CD14+ Monocytes
4  AAACCGTGCTTCCG-1   CD79A    0.0  donor_1     control  CD14+ Monocytes
(6351, 6)
```

Seaborn consumes `plot_data` directly:

```python
import seaborn as sns

sns.catplot(
    data=plot_data,
    x="louvain",
    y="value",
    hue="condition",
    col="feature",
    kind="strip",
    sharey=False,
)
```

The plotting dependency stays outside annplyr, and the table contains exactly
the 6,351 matrix values requested. The finite budget documents that expectation
and fails an oversized request before the first matrix read.

## Keep domain-specific steps composable

`pipe()` is the escape hatch for project-specific functions — take AnnData, use
ordinary Python or another scverse package, return whatever the next step needs:

```{testcode}
def donor_counts(data):
    return data.ap.count("donor_id", sort=True)


print(analysis.ap.pipe(donor_counts))
```

```{testoutput}
  donor_id    n
0  donor_1  711
1  donor_3  709
2  donor_2  697
```

`analysis` is still a standard AnnData object, ready for Scanpy, CellRank,
scvi-tools, serialization, or another `.ap` pipeline. annplyr owns the wrangling
semantics; the surrounding workflow keeps ownership of normalization,
modelling, visualization, and biological interpretation.

## The whole workflow

The central steps read as one pipeline:

```python
analysis = (
    adata.ap.left_join(sample_sheet, by="donor_id", relationship="many-to-one")
    .ap.filter(
        obs=(
            ap.col("n_genes") >= 500,
            ap.col("percent_mito") < 0.03,
        )
    )
    .ap.mutate(
        raw={"b_score": (ap.col("MS4A1") + ap.col("CD79A")) / 2},
        max_matrix_values=2 * adata.n_obs,
    )
    .ap.select(
        obs=["donor_id", "condition", "louvain", "b_score"],
        x=["MS4A1", "CD79A", "NKG7"],
    )
)
```

Each call stays an ordinary, inspectable operation, but the pipeline states the
analysis more directly than a sequence of boolean masks, `.obs` assignments,
pandas merges, and hand-synchronized AnnData slices.
