# Single-Cell Utilities

`annplyr` includes small AnnData utility helpers for metadata workflows. These
helpers are intentionally narrow; they are not a biological QC or annotation
toolkit.

## Sample Metadata

```python
meta = ap.sample_meta(adata, sample="sample_id", include=["batch", "condition"])
summary = ap.sample_summary(adata, sample="sample_id")
adata = ap.add_sample_meta(adata, sample="sample_id", meta=sample_table)
```

## Feature Presence

```python
presence = ap.feature_present(adata, ["MS4A1", "CD79A", "missing_gene"])
presence.to_frame()
```

## Safe Names

```python
adata = ap.rename_obs_names(adata, lambda name: f"cell_{name}")
adata = ap.add_name_prefix(adata, "sample1", axis="obs")
duplicates = ap.name_duplicates(adata, axis="obs")
```

## Palettes

Palettes are stored using Scanpy-compatible `uns` keys.

```python
adata = ap.store_palette(adata, "cell_type", ["#1f77b4", "#ff7f0e"])
palette = ap.get_palette(adata, "cell_type")
```

Out of scope for core `annplyr`: species gene registries, mitochondrial or
ribosomal scoring, cell-cycle scoring, QC plotting wrappers, marker discovery,
and cluster annotation.
