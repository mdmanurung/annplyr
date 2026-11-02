# Safety Notes

## AnnData Alignment

AnnData has coordinated containers. If observations or features are subset or reordered, the matching axis must stay aligned in `X`, layers, `obsm`, `varm`, `obsp`, and `varp`. Prefer annplyr accessor verbs or AnnData-native slicing instead of manual pandas filtering plus reassignment.

Every verb returns a new object, so assign the result. Subsetting verbs return an independent AnnData; a backed input returns the view. Same-shape verbs copy `obs`, `var` and `uns` and share `X`, layers and the aligned matrices with the input.

## Matrix Materialization

Sparse matrices and backed arrays can be large. Predicates over selected matrix columns are fine; whole-matrix pandas exports should be explicit. For long exports, name the features; a whole `X`, `raw` or layer needs `allow_all_features=True`, which the user must intend. Use cumulative `max_matrix_values=` on matrix-reading verbs; annplyr plans every source and rejects over-budget requests before the first read.

Canonical scalar matrix summaries automatically use deterministic feature or row chunks across dense, CSR/CSC, and backed sources. There is no public chunk-size flag. Exact median and distinct-count reducers can retain one reduction vector's state; use a finite cumulative budget for opaque or cross-row custom summary expressions that require the eager compatibility path.

## Mutation

`mutate()` and `transmute()` write metadata columns in `obs` or `var`. They may read from `X`, layers, `raw`, `obsm`, and `varm`, but they should not modify those matrices.

## Controlled Extraction

Use `pull()` for inspection of non-core AnnData containers such as `raw`, `obsp`, `varp`, or tabular `uns`. Pairwise matrices are extraction-only; do not try to join or mutate AnnData axes from them.

## Backed Objects

A subsetting verb on a backed input returns the view, so a large file can be subset without loading it. Backed objects cannot support metadata-changing verbs; load into memory first when that is intentional:

```python
adata = adata.to_memory()
adata = adata.mutate(new_column=obs.old_column)
```

## Good Failure Modes

Prefer annplyr's typed errors over broad `ValueError` paths in package code:

- `SelectionError` for invalid selectors, and when `relocate`'s `before`/`after` anchor is among the columns being moved.
- `UnknownColumnError` or `UnknownSourceError` for unresolved expression sources.
- `DuplicateNameError` and `NameRepairError` for schema collisions.
- `IncompatibleAxisError` and `SizeMismatchError` for axis or length mismatches.
- `JoinRelationshipError` for unsafe join cardinality.

## Parsing and Nesting Helpers

`ap.tidyr.separate(df, col, into=..., sep=..., ...)` always treats `sep` as a regular expression. NA values in the source column yield `[pd.NA] * len(into)`; the string `"None"` is never split.

`ap.tidyr.unnest(df, col)` captures inner column names from the first non-empty nested frame, so an all-empty result still has the correct schema (the expected columns are present but with zero rows).
