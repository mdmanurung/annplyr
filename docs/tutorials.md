# Tutorials

The tutorials use deterministic, synthetic single-cell data so every notebook
is redistributable, quick to run, and focused on workflow structure. The
biological labels are realistic, but the values are illustrative rather than a
reference dataset.

```{toctree}
:maxdepth: 1

vignettes/cohort_workflow
notebooks/getting_started
notebooks/all_of_annplyr
notebooks/plot_ready_tables
```

## Choose a tutorial

1. {doc}`vignettes/cohort_workflow` tells the full story from an annotated
   AnnData object to sample-aware summaries and plot-ready marker data.
2. {doc}`notebooks/getting_started` is a self-contained notebook version of a
   compact cohort analysis, including sample-sheet joins and grouped ranking.
3. {doc}`notebooks/all_of_annplyr` is a cookbook for selectors, expressions,
   core verbs, grouping, extraction, rectangling, and `pipe()`.
4. {doc}`notebooks/plot_ready_tables` focuses on tidy expression tables,
   sample-level summaries, and model-ready wide tables.

For a specific method rather than a narrative workflow, use the {doc}`api`.
