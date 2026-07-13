# Joins

Joins in `annplyr` are AnnData-safe metadata joins. They can enrich or subset
`obs` or `var`, but they cannot silently create new cells, new features, or
duplicated axis records.

## Observation Metadata

```python
joined = adata.ap.left_join(sample_table, by="sample_id", axis="obs")
```

## Feature Metadata

```python
annotated = adata.ap.left_join(gene_table, by="gene_id", axis="var")
```

## Filtering Joins

```python
subset = adata.ap.semi_join(allowed_samples, by="sample_id", axis="obs")
removed = adata.ap.anti_join(excluded_genes, by="gene_id", axis="var")
```

## Relationship Checks

Use relationship arguments when cardinality matters:

```python
adata.ap.left_join(
    sample_table,
    by="sample_id",
    axis="obs",
    relationship="many-to-one",
)
```

If a join would duplicate or add AnnData axis records, `annplyr` raises
`JoinRelationshipError`.
