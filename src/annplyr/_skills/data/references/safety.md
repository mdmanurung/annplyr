# Safety Notes

## AnnData Alignment

AnnData has coordinated containers. If observations or features are subset or reordered, the matching axis must stay aligned in `X`, layers, `obsm`, `varm`, `obsp`, and `varp`. Prefer annplyr accessor verbs or AnnData-native slicing instead of manual pandas filtering plus reassignment.

## Matrix Materialization

Sparse matrices and backed arrays can be large. Predicates over selected matrix columns are fine; whole-matrix pandas exports should be explicit. For long exports, pass selected features unless the user intentionally accepts all-feature materialization. Use `max_matrix_values=` on `to_df()`, `to_tidy()`, `pivot_longer()`, or `as_frame()` when a workflow needs a hard materialization budget.

## Mutation

`mutate()` and `transmute()` write metadata columns in `obs` or `var`. They may read from `x`, layers, `raw`, `obsm`, and `varm`, but they should not modify those matrices.

## Controlled Extraction

Use `as_frame()` for inspection of non-core AnnData containers such as `raw`, `obsp`, `varp`, or tabular `uns`. Pairwise matrices are extraction-only; do not try to join or mutate AnnData axes from them.

## Backed Objects

Backed AnnData objects cannot safely support metadata mutation through annplyr. Load into memory first when mutation is intentional:

```python
adata = adata.to_memory()
adata = adata.ap.mutate(obs={"new_column": ap.col("old_column")})
```

## Good Failure Modes

Prefer annplyr's typed errors over broad `ValueError` paths in package code:

- `SelectionError` for invalid selectors.
- `UnknownColumnError` or `UnknownSourceError` for unresolved expression sources.
- `DuplicateNameError` and `NameRepairError` for schema collisions.
- `IncompatibleAxisError` and `SizeMismatchError` for axis or length mismatches.
- `JoinRelationshipError` for unsafe join cardinality.
