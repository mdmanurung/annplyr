# annplyr v0.3 public API contract

Status: frozen for v0.3.0 implementation on 2026-08-08.

This inventory is normative. The first column is a machine-readable identifier;
`tests/test_api_contract_v03.py` requires one row for every accessor method,
grouped method/introspection operation, and exported public symbol. No cell may
be empty or contain a placeholder.

Exact signatures, overloads, return annotations, and the public
`annplyr.typing` aliases are generated in {doc}`public-typing-contract`; both
inventories must remain current.

## Contract vocabulary

- **shape-copy**: the result is independent mutable storage across `X`,
  layers, raw, `obsm`, `varm`, `obsp`, and `varp`. There is no `copy` option.
  A backed input returns the backed view, so a backed `.h5ad` can be subset
  without loading it.
- **same-shape**: the result is a new object that shares `X`, layers, raw,
  `obsm`, `varm`, `obsp`, and `varp` with the input, with copied `obs`, `var`,
  and `uns`. There is no `inplace` option. Backed input fails before writing.
- **independent**: the result is independent and in memory.
- **positional**: axis identity is carried only by zero-based integer positions;
  duplicate axis names are never used for reconstruction. Survivors retain
  input order unless the verb explicitly sorts.
- **projected budget**: `max_matrix_values=None` is unbounded. A nonnegative
  value is charged by cumulative logical projected matrix cells across the
  whole call before any read. Opaque expressions charge their full source.
  Invalid, negative, or cumulative over-budget calls perform zero reads.
- **container reference**: an expression names the container it reads from.
  A bare name is `obs` metadata. A matrix always needs a qualifier, because no
  default is safe: `x:` for `X`, `raw:` for `.raw`, a layer name such as
  `counts:`, an `obsm`/`varm` key such as `X_pca:`, or `var:`/`varm:` for the
  feature axis. `obs:` and `var:` are the explicit forms of the bare namespace.
  A column name that legitimately contains a colon, such as `ENSG001:1`, is
  never reinterpreted, because only a reserved or known-container prefix
  qualifies.
- **one container per expression**: an expression whose references span two
  containers is rejected with `IncompatibleAxisError`, because each container is
  evaluated against its own frame. `mutate` writes whichever container the
  caller names through `into=`, so assignments may read from any container. An
  expression naming no column, such as a literal or a raw Narwhals expression,
  is evaluated against the default `obs` frame unless the verb's `source=` names
  a container, in which case it is charged for that container in full.
- **grouping key**: a grouping key names metadata on the reduced axis, and what
  it must say depends on whether the verb already knows that axis. `group_by`
  takes an explicit `axis=`, so a key is a bare name on that axis and a
  qualifier is redundant. `summarize` and `count` learn the axis from the
  aggregates, so a bare name is `obs` and a `var:` name is a feature, making the
  qualifier required to disambiguate a feature key. A key naming the other axis
  is rejected with `IncompatibleAxisError` rather than resolved, because
  resolving `var:kind` against an `obs` frame would silently group whichever
  `obs` column shares the bare name.
- **schema-only**: selectors resolve from names and dtypes without reading
  matrix values. Matrix `where()` receives a zero-length typed Series.
- **grouped result**: an AnnData-returning grouped call returns a new
  `GroupedAnnData` with an updated `GroupSpec`. `ungroup()` returns the exact underlying object, whose
  class is upgraded in place as described under **chainable result**.
- **group plan**: one positional plan per public call; observed nonempty groups
  include NA keys and appear in first-seen order. `.rows` stores zero-based
  positions and key columns retain extension/categorical dtypes.
- **chainable result**: every AnnData-returning accessor verb returns a
  `TidyAnnData`, which is an `AnnData` subclass carrying every verb as a
  method, so `isinstance(x, AnnData)` still holds and a chain needs `.ap` only
  once. The upgrade is an in-place class change, so object identity is
  preserved.
  When the input is a user-defined `AnnData` subclass, the upgrade composes a
  class with it, preserving that subclass's methods, slots, and `isinstance`
  relationships. Verbs that derive a new object (`filter`,
  `select`, `mutate`, and friends) return a plain `TidyAnnData`, matching upstream
  `AnnData`, whose `copy` and `__getitem__` also return a plain `AnnData` for
  subclasses. `TidyAnnData.to_anndata()` returns an independent `AnnData`
  without the verbs, restoring a user subclass the object still carried.

Typed failures use the package hierarchy: `SelectionError` for invalid
selectors, budgets, or scalar parameters; `UnknownSourceError` and
`UnknownColumnError` for missing inputs; `IncompatibleAxisError` for illegal
axis/source combinations; `SizeMismatchError` for invalid result sizes;
`DuplicateNameError` or `NameRepairError` for ambiguous output schemas;
`JoinRelationshipError` for join cardinality/right-only/unmatched violations;
and `AnnplyrError` for backed writes or invalid pipe targets.

Post-v0.3 readiness addendum: canonical scalar matrix summaries retain this
public contract while executing through deterministic feature/row chunks.
Their complete cumulative budget is still validated before the first adapter
read; result values, dtypes, key order, and positional identity match the eager
contract. Opaque and cross-row custom expressions keep their documented
full-source conservative charge.

## AnnData accessor inventory

| ID | Axes | Accepted sources | Return | Grouping | Ordering | Ownership | Sparse/backed | Budget | Typed failures |
|---|---|---|---|---|---|---|---|---|---|
| `accessor.filter` | obs and var | one predicate per argument, each naming its own container | AnnData | ungrouped | positional survivor order | shape-copy | preserve sparse format; backed input returns the view | projected budget | selection, source, column, axis, budget, size |
| `accessor.select` | obs metadata, var metadata, X columns | obs, var, X schema selectors | AnnData | ungrouped | requested column order; positional var order for X | shape-copy | preserve sparse format; backed input returns the view | schema-only; no charge | selection, source, column, duplicate name |
| `accessor.rename` | obs metadata, var metadata, X columns | explicit old-to-new mappings | AnnData | ungrouped | stable except requested names | same-shape | preserve sparse format; backed write rejected | schema-only; no charge | column, duplicate name, backed write |
| `accessor.rename_with` | obs metadata, var metadata, X columns | callable plus schema selectors | AnnData | ungrouped | stable except requested names | same-shape | preserve sparse format; backed write rejected | schema-only; no charge | selection, column, duplicate name, backed write |
| `accessor.relocate` | obs metadata, var metadata, X columns | schema selectors and before/after anchor | AnnData | ungrouped | requested columns relocated; X relocation reorders var positionally | same-shape | preserve sparse format; backed write rejected | schema-only; no charge | selection, column, axis, backed write |
| `accessor.distinct` | selected axis | obs, var, X expressions | AnnData | ungrouped | positional first occurrence | shape-copy | preserve sparse format; backed input returns the view | projected budget | selection, source, column, axis, budget |
| `accessor.left_join` | selected axis | pandas DataFrame or mapping metadata | AnnData | ungrouped | left positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | join relationship, column, name repair, axis |
| `accessor.inner_join` | selected axis | pandas DataFrame or mapping metadata | AnnData | ungrouped | surviving left positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | join relationship, column, name repair, axis |
| `accessor.right_join` | selected axis | pandas DataFrame or mapping metadata | AnnData | ungrouped | right-key order then left order within key | shape-copy | preserve sparse format; reject right-only records | no matrix reads | join relationship, right-only, column, name repair, axis |
| `accessor.full_join` | selected axis | pandas DataFrame or mapping metadata | AnnData | ungrouped | left positional order | shape-copy | preserve sparse format; reject right-only records | no matrix reads | join relationship, right-only, column, name repair, axis |
| `accessor.semi_join` | selected axis | pandas DataFrame or mapping metadata | AnnData | ungrouped | surviving left positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | column, axis, NA mode |
| `accessor.anti_join` | selected axis | pandas DataFrame or mapping metadata | AnnData | ungrouped | surviving left positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | column, axis, NA mode |
| `accessor.arrange` | obs and var | one sort key per argument, each naming its own container | AnnData | ungrouped | stable explicit sort, positional ties | shape-copy | preserve sparse format; backed input returns the view | projected budget | selection, source, column, axis, budget |
| `accessor.slice` | selected axis | integer positions and slices | AnnData | ungrouped | requested positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | selection, axis |
| `accessor.slice_head` | selected axis | nonnegative n | AnnData | ungrouped | first n positions | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | selection, axis |
| `accessor.slice_tail` | selected axis | nonnegative n | AnnData | ungrouped | last n positions | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | selection, axis |
| `accessor.slice_min` | selected axis | metadata expression and n | AnnData | ungrouped | stable ascending expression order | shape-copy | preserve sparse format; backed input returns the view | no matrix source reads | selection, column, axis |
| `accessor.slice_max` | selected axis | metadata expression and n | AnnData | ungrouped | stable descending expression order | shape-copy | preserve sparse format; backed input returns the view | no matrix source reads | selection, column, axis |
| `accessor.slice_sample` | selected axis | n or prop, replace, random state | AnnData | ungrouped | RNG selection order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | selection, axis |
| `accessor.mutate` | obs and var metadata | one assignment dict; each value names its own container, `into=` picks the written axis | AnnData | ungrouped | axis order stable; assignments strictly sequential | same-shape | preserve sparse format; backed write rejected | projected budget | selection, source, column, size, budget, backed write |
| `accessor.transmute` | obs and var metadata | one assignment dict; each value names its own container, `into=` picks the written axis | AnnData | ungrouped | requested assignment order; sequential dependencies | independent | preserve sparse format; backed input yields in-memory result | projected budget | selection, source, column, size, budget |
| `accessor.group_by` | exactly one of obs or var, given as `axis=` | existing metadata column names or tidy schema selectors, one per positional key | GroupedAnnData or exact AnnData when no axis supplied | persistent GroupSpec | group plan first-seen order | wrapper shares exact underlying object | no matrix, virtual, computed, or empty grouping | no matrix reads | selection, column, incompatible axis |
| `accessor.summarize` | obs and var evaluation axes | one aggregate dict; each aggregate names its own container | pandas DataFrame | optional ephemeral by selector | first-seen key order | result table independent | sparse/backed projected reads supported | projected budget | selection, source, column, size, budget |
| `accessor.summarise` | obs and var evaluation axes | same as summarize | pandas DataFrame | optional ephemeral by selector | first-seen key order | alias of summarize | sparse/backed projected reads supported | projected budget | same as summarize |
| `accessor.count` | selected metadata axis | metadata grouping selectors and optional weight | pandas DataFrame | ephemeral by selector | first-seen or descending count when sort | result table independent | metadata only; backed safe | no matrix reads | selection, column, axis |
| `accessor.add_count` | selected metadata axis | metadata grouping selectors and optional weight | AnnData | ungrouped | stable unless sort requests axis reorder | same-shape | backed write rejected; sparse format preserved | no matrix reads | selection, column, axis, size, backed write |
| `accessor.pull` | source-dependent | exactly one column from obs, var, X/layer, raw, one obsm/varm/obsp/varp key, or one uns key | pandas Series | ungrouped | source positional order | independent pandas result | dense, sparse, DataFrame, and backed projections | projected budget | selection, source, column, budget |
| `accessor.glimpse` | obs rows and container names | row count | None | ungrouped | obs row order | no result object | prints only; no AnnData change | no matrix source reads | none |
| `accessor.to_wide` | obs rows | obs, X/layer, raw, obsm, obsp | pandas DataFrame | ungrouped | obs positional order and requested columns | independent pandas result | wide sparse uses pandas sparse columns; backed projection | projected budget | selection, source, column, name repair, budget |
| `accessor.join_features` | obs rows by features | feature references from X, raw or a layer; shape long or wide | pandas DataFrame (long: obs columns plus `.cell`, `.feature`, `.abundance_<source>`) or AnnData (wide) | ungrouped | obs-major then requested feature order | independent long table; wide is a new AnnData | only selected dense/sparse/backed values read | projected budget | selection, source, column, axis, name repair, duplicate name, budget |
| `accessor.pivot_longer` | obs rows by features | obs, X/layer, raw | pandas DataFrame | ungrouped | obs-major then selected feature order | independent long table | only selected dense/sparse/backed values read | projected budget | selection, source, column, name repair, budget |
| `accessor.aggregate_cells` | obs rows | obs key selectors, required source matrix, fun sum or mean | AnnData | one row per observed key combination | first-seen key order | new AnnData with keys and `n_cells` in obs | sparse one-hot product; reads the whole source | no projected budget | selection, source, column, axis |
| `accessor.nest` | obs rows by group | obs key selectors, output name | pandas DataFrame of keys and AnnData subsets | keys emitted | first-seen key order | each subset is an independent AnnData | subset by position | no matrix evaluation reads | selection, column, name |
| `accessor.pipe` | generic | callable or callable/keyword tuple | callable-defined | passes exact AnnData | callable-defined | callable-defined | callable-defined | callable-defined | duplicate pipe keyword, invalid callable |

## Persistent grouped inventory

All AnnData-returning rows below use **grouped result** semantics. Group-local
shape verbs concatenate groups in first-seen order and retain positional order
within each group unless explicitly sorted. Joins execute exactly once globally.

| ID | Axes | Accepted sources | Return | Grouping | Ordering | Ownership | Sparse/backed | Budget | Typed failures |
|---|---|---|---|---|---|---|---|---|---|
| `grouped.__iter__` | GroupSpec axis | stored metadata keys | iterator of key dict and chainable AnnData view | persistent | group plan order | group views, wrapper unchanged | sparse/backed views allowed | no reads beyond metadata | invalid stored grouping |
| `grouped.group_vars` | GroupSpec axis | stored metadata keys | list of strings | introspection | GroupSpec column order | independent list | metadata only | no matrix reads | invalid stored grouping |
| `grouped.group_keys` | GroupSpec axis | stored metadata keys | pandas DataFrame | introspection | group plan order | independent table; key dtypes preserved | metadata only | no matrix reads | invalid stored grouping |
| `grouped.group_data` | GroupSpec axis | stored metadata keys | pandas DataFrame with `.rows` | introspection | group plan order and zero-based positional row lists | independent table; key dtypes preserved | metadata only | no matrix reads | invalid stored grouping |
| `grouped.ungroup` | GroupSpec axis | none | exact underlying AnnData, upgraded in place to a chainable result | removed | unchanged | identity return; class upgraded in place | unchanged | no reads | none |
| `grouped.filter` | both, group-local on GroupSpec axis | one predicate per argument, each naming its own container | GroupedAnnData | preserved and regrouped from result | group plan then positional survivors | grouped shape-copy | sparse preserved; backed input returns the view | projected budget | accessor filter failures plus missing keys |
| `grouped.select` | obs metadata, var metadata, X columns | accessor select sources; keys auto-retained | GroupedAnnData | keys retained, omitted keys placed first | requested order with deterministic key insertion | grouped shape-copy | sparse preserved; backed input returns the view | schema-only | accessor select failures plus missing keys |
| `grouped.rename` | obs metadata, var metadata, X columns | explicit mappings | GroupedAnnData | GroupSpec names updated | stable | grouped same-shape | backed write rejected; sparse preserved | schema-only | accessor rename failures plus ambiguous key output |
| `grouped.rename_with` | obs metadata, var metadata, X columns | callable plus selectors | GroupedAnnData | GroupSpec names updated | stable | grouped same-shape | backed write rejected; sparse preserved | schema-only | accessor rename failures plus ambiguous key output |
| `grouped.relocate` | obs metadata, var metadata, X columns | accessor relocate sources | GroupedAnnData | preserved | requested relocation | grouped same-shape | backed write rejected; sparse preserved | schema-only | accessor relocate failures plus missing keys |
| `grouped.transmute` | obs and var metadata | one assignment dict; keys auto-retained | GroupedAnnData | keys retained and final values regroup future calls | group plan evaluation then assignment order | grouped independent | sparse preserved; backed result in memory | projected budget | accessor transmute failures plus missing keys |
| `grouped.left_join` | selected join axis | pandas DataFrame or mapping metadata | GroupedAnnData | retained; left collision suffix updates GroupSpec | global left order | grouped shape-copy | sparse preserved; backed input returns the view | no matrix reads | accessor join failures plus lost/ambiguous grouping key |
| `grouped.inner_join` | selected join axis | pandas DataFrame or mapping metadata | GroupedAnnData | retained; left collision suffix updates GroupSpec | global surviving left order | grouped shape-copy | sparse preserved; backed input returns the view | no matrix reads | accessor join failures plus lost/ambiguous grouping key |
| `grouped.right_join` | selected join axis | pandas DataFrame or mapping metadata | GroupedAnnData | retained; left collision suffix updates GroupSpec | global right-key then left-within-key order | grouped shape-copy | sparse preserved; right-only rejected | no matrix reads | accessor join failures plus lost/ambiguous grouping key |
| `grouped.full_join` | selected join axis | pandas DataFrame or mapping metadata | GroupedAnnData | retained; left collision suffix updates GroupSpec | global left order | grouped shape-copy | sparse preserved; right-only rejected | no matrix reads | accessor join failures plus lost/ambiguous grouping key |
| `grouped.semi_join` | selected join axis | pandas DataFrame or mapping metadata | GroupedAnnData | retained | global surviving left order | grouped shape-copy | sparse preserved; backed input returns the view | no matrix reads | accessor join failures plus lost grouping key |
| `grouped.anti_join` | selected join axis | pandas DataFrame or mapping metadata | GroupedAnnData | retained | global surviving left order | grouped shape-copy | sparse preserved; backed input returns the view | no matrix reads | accessor join failures plus lost grouping key |
| `grouped.mutate` | both, group-local on GroupSpec axis | one assignment dict; keys auto-retained | GroupedAnnData | call uses old plan; final key values regroup future calls | original axis order; sequential within group | grouped same-shape | backed write rejected; sparse preserved | projected budget | accessor mutate failures plus dtype reconstruction |
| `grouped.add_count` | GroupSpec axis | optional metadata weight | GroupedAnnData | preserved and regrouped after call | stable unless sort reorders groups | grouped same-shape | backed write rejected; sparse preserved | no matrix reads | selection, size, backed write |
| `grouped.arrange` | both, group-local on GroupSpec axis | one sort key per argument, each naming its own container | GroupedAnnData | preserved | group plan then stable local sort | grouped shape-copy | sparse preserved; backed input returns the view | projected budget | accessor arrange failures plus missing keys |
| `grouped.distinct` | GroupSpec axis | matching-axis metadata and X where compatible | GroupedAnnData | preserved | group plan then local first occurrence | grouped shape-copy | sparse preserved; backed input returns the view | projected budget | accessor distinct failures plus incompatible axis |
| `grouped.slice` | GroupSpec axis | integer local positions and slices | GroupedAnnData | preserved | group plan then requested local positions | grouped shape-copy | sparse preserved; backed input returns the view | no matrix reads | selection |
| `grouped.slice_head` | GroupSpec axis | nonnegative n | GroupedAnnData | preserved | group plan then first n local positions | grouped shape-copy | sparse preserved; backed input returns the view | no matrix reads | selection |
| `grouped.slice_tail` | GroupSpec axis | nonnegative n | GroupedAnnData | preserved | group plan then last n local positions | grouped shape-copy | sparse preserved; backed input returns the view | no matrix reads | selection |
| `grouped.slice_min` | GroupSpec axis | local metadata expression and n | GroupedAnnData | preserved | group plan then stable ascending local sort | grouped shape-copy | sparse preserved; backed input returns the view | no matrix source reads | selection, column |
| `grouped.slice_max` | GroupSpec axis | local metadata expression and n | GroupedAnnData | preserved | group plan then stable descending local sort | grouped shape-copy | sparse preserved; backed input returns the view | no matrix source reads | selection, column |
| `grouped.slice_sample` | GroupSpec axis | n or prop, replace, random state | GroupedAnnData | preserved | group plan then RNG local selection | grouped shape-copy | sparse preserved; backed input returns the view | no matrix reads | selection |
| `grouped.summarize` | both, group-local on GroupSpec axis | one aggregate dict, grouped by the stored keys | pandas DataFrame | keys emitted | group plan order | independent table; key dtypes preserved | dense, sparse, DataFrame, backed projection | projected budget | accessor summarize failures |
| `grouped.summarise` | both, group-local on GroupSpec axis | same as grouped summarize | pandas DataFrame | keys emitted | group plan order | alias of summarize | dense, sparse, DataFrame, backed projection | projected budget | same as grouped summarize |
| `grouped.count` | GroupSpec axis | optional metadata weight | pandas DataFrame | keys emitted | group plan or descending count | independent table; key dtypes preserved | metadata only | no matrix reads | selection, column |
| `grouped.aggregate_cells` | GroupSpec axis | required source matrix, fun sum or mean | AnnData | group keys become the rows | group plan order | new AnnData with keys and `n_cells` in obs | sparse one-hot product; reads the whole source | no projected budget | selection, source, column |
| `grouped.pipe` | generic | callable or callable/keyword tuple | callable-defined | wrapper is the pipe target | callable-defined | callable-defined | callable-defined | callable-defined | duplicate pipe keyword, invalid callable |

## Exported symbol inventory

Expression helpers return `AnnplyrExpr` unless a row says selector/helper.
`AnnplyrExpr` exposes `to_narwhals()` and proxies public operators, methods, and
namespaces while conservatively preserving exact/opaque dependencies, output
width, and scalar/row/unknown cardinality. Raw Narwhals expressions remain
accepted as isolated full-source opaque inputs.

| ID | Axes | Accepted sources | Return | Grouping | Ordering | Ownership | Sparse/backed | Budget | Typed failures |
|---|---|---|---|---|---|---|---|---|---|
| `export.AnnplyrAccessor` | both | one AnnData instance | registered accessor instance | ungrouped entry point | delegates to each verb | cached by AnnData namespace registration | delegates to each verb | delegates to each verb | constructor signature and delegated typed failures |
| `export.AnnplyrError` | n/a | message | exception instance | n/a | n/a | independent | n/a | n/a | base typed package error |
| `export.AnnplyrExpr` | schema-dependent | wrapped Narwhals expression and dependency metadata | immutable expression wrapper | caller-defined | expression-defined | independent helper | adapter chosen by caller | caller projected budget | expression, selection, and size failures |
| `export.DuplicateNameError` | n/a | message | exception instance | n/a | n/a | independent | n/a | n/a | duplicate output names |
| `export.FeaturePresence` | n/a | found and missing collections | report object with mutable list fields | n/a | input order | independent | n/a | n/a | constructor validation |
| `export.GroupedAnnData` | obs or var | AnnData and one resolved grouping selector | persistent grouped wrapper | retains positional GroupSpec | stable group plan order | wrapper owns the referenced AnnData result | delegates to grouped verbs | delegates to grouped verbs | selection, axis, join, size, and budget failures |
| `export.IncompatibleAxisError` | n/a | message | exception instance | n/a | n/a | independent | n/a | n/a | incompatible axes |
| `export.JoinRelationshipError` | n/a | message | exception instance | n/a | n/a | independent | n/a | n/a | join cardinality/unmatched/right-only |
| `export.NameRepairError` | n/a | message | exception instance | n/a | n/a | independent | n/a | n/a | unsafe name repair |
| `export.SelectionError` | n/a | message | exception instance | n/a | n/a | independent | n/a | n/a | invalid selection or budget |
| `export.SizeMismatchError` | n/a | message | exception instance | n/a | n/a | independent | n/a | n/a | invalid output cardinality |
| `export.TidyAnnData` | both | an existing AnnData object | AnnData subclass carrying every verb as a method | grouping stays on GroupedAnnData | delegates to each verb | wrap preserves object identity | delegates to each verb | delegates to each verb | constructor accepts an AnnData; to_anndata drops the verbs and restores any carried user subclass |
| `export.UnknownColumnError` | n/a | message | exception instance | n/a | n/a | independent | n/a | n/a | missing column |
| `export.UnknownSourceError` | n/a | message | exception instance | n/a | n/a | independent | n/a | n/a | missing source |
| `export.across` | schema-dependent | tidy selector and expression functions | assignment expander carrying wrapped expressions | caller-defined | selector order | immutable helper | schema-only resolution; adapter reads by caller | caller projected budget | selection, duplicate output name |
| `export.add_name_prefix` | obs or var names | AnnData, prefix, axis | AnnData | ungrouped | stable | same-shape | backed write rejected; sparse unchanged | no matrix reads | axis, duplicate name, backed write |
| `export.all_of` | schema-dependent | required names | tidy selector | caller-defined | requested order | immutable helper | schema-only | no direct read | unknown column, selection |
| `export.any_of` | schema-dependent | optional names | tidy selector | caller-defined | requested surviving order | immutable helper | schema-only | no direct read | selection |
| `export.average_rank` | row expression | string or wrapped/raw expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.between` | row expression | expression and bounds | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.case_match` | row expression | expression, cases, default | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.case_when` | row expression | predicate/value cases and default | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.coalesce` | row expression | expressions or literals | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.col` | schema-dependent | one or more column names | AnnplyrExpr with exact dependencies | caller-defined | requested expression order | immutable expression | adapter chosen by caller | caller projected budget | unknown column |
| `export.contains` | schema-dependent | substring and regex flag | tidy selector | caller-defined | schema order | immutable helper | schema-only | no direct read | invalid pattern |
| `export.cum_max` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.cum_min` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.cum_prod` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.cum_sum` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.cumall` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.cumany` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.cume_dist` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.cummean` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.dense_rank` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.desc` | row expression | string or expression | descending-order helper carrying AnnplyrExpr | caller-defined | descending marker | immutable helper | adapter chosen by caller | caller projected budget | expression |
| `export.ends_with` | schema-dependent | suffix | tidy selector | caller-defined | schema order | immutable helper | schema-only | no direct read | selection |
| `export.everything` | schema-dependent | none | tidy selector | caller-defined | schema order | immutable helper | schema-only | no direct read | none |
| `export.feature_present` | var | AnnData and feature names | FeaturePresence | none | requested feature order | independent report | metadata/index only; backed safe | no matrix reads | selection |
| `export.first` | aggregate expression | string or expression | AnnplyrExpr | caller-defined | first row | immutable expression | adapter chosen by caller | caller projected budget | empty/value expression |
| `export.get_palette` | obs | AnnData and categorical column/key | color mapping | none | category order | independent mapping | metadata/uns only; backed safe | no matrix reads | unknown column, palette mismatch |
| `export.glimpse` | obs rows and container names | AnnData and row count | None | ungrouped | obs row order | no result object | prints only; no AnnData change | no matrix source reads | none |
| `export.if_all` | row expression | selector and predicate | wrapped expression expander | caller-defined | row order | immutable helper | schema projection by caller | caller projected budget | selection, expression |
| `export.if_any` | row expression | selector and predicate | wrapped expression expander | caller-defined | row order | immutable helper | schema projection by caller | caller projected budget | selection, expression |
| `export.if_else` | row expression | condition and two values | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.join_features` | obs rows by features | AnnData, feature references, shape; adata passed as the first argument | pandas DataFrame or AnnData | ungrouped | obs-major then requested feature order | independent long table; wide is a new AnnData | only selected dense/sparse/backed values read | projected budget | selection, source, column, axis, name repair, duplicate name, budget |
| `export.is_na` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression |
| `export.lag` | row expression | expression, offset, default | AnnplyrExpr | caller-defined | input row order shifted backward | immutable expression | adapter chosen by caller | caller projected budget | offset, expression, size |
| `export.last` | aggregate expression | string or expression | AnnplyrExpr | caller-defined | last row | immutable expression | adapter chosen by caller | caller projected budget | empty/value expression |
| `export.last_col` | schema-dependent | nonnegative offset | tidy selector | caller-defined | selected schema position | immutable helper | schema-only | no direct read | selection |
| `export.layer` | container name | layer, obsm key, or varm key name | Source namespace | caller-defined | caller-defined | immutable helper | no source read itself | no direct read | none |
| `export.lead` | row expression | expression, offset, default | AnnplyrExpr | caller-defined | input row order shifted forward | immutable expression | adapter chosen by caller | caller projected budget | offset, expression, size |
| `export.lit` | row/scalar expression | literal value | AnnplyrExpr with no dependencies | caller-defined | caller-defined | immutable expression | no source read itself | caller projected budget | unsupported literal |
| `export.Source` | container name | prefix of one container | namespace building column expressions | caller-defined | caller-defined | immutable helper | no source read itself | no direct read | TypeError on call |
| `export.matches` | schema-dependent | regex | tidy selector | caller-defined | schema order | immutable helper | schema-only | no direct read | invalid regex |
| `export.max` | aggregate expression | string or expression | AnnplyrExpr | caller-defined | scalar | immutable expression | adapter chosen by caller | caller projected budget | expression |
| `export.max_rank` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.mean` | aggregate expression | string or expression | AnnplyrExpr | caller-defined | scalar | immutable expression | adapter chosen by caller | caller projected budget | expression |
| `export.median` | aggregate expression | string or expression | AnnplyrExpr | caller-defined | scalar | immutable expression | adapter chosen by caller | caller projected budget | expression |
| `export.min` | aggregate expression | string or expression | AnnplyrExpr | caller-defined | scalar | immutable expression | adapter chosen by caller | caller projected budget | expression |
| `export.min_rank` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.n` | aggregate expression | none | AnnplyrExpr with no column dependencies | caller-defined | scalar | immutable expression | no source column itself | caller projected budget | expression context |
| `export.n_distinct` | aggregate expression | string or expression | AnnplyrExpr | caller-defined | scalar | immutable expression | adapter chosen by caller | caller projected budget | expression |
| `export.na_if` | row expression | expression and comparison value | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, dtype |
| `export.name_duplicates` | obs or var names | AnnData and axis | pandas DataFrame | none | first duplicate position order | independent table | index only; backed safe | no matrix reads | axis |
| `export.near` | row expression | expression, other value, tolerance | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | tolerance, expression |
| `export.nest` | obs rows by group | AnnData and obs key selectors, output name; adata passed as the first argument | pandas DataFrame of keys and AnnData subsets | keys emitted | first-seen key order | each subset is an independent AnnData | subset by position; backed input is materialized per group | no matrix evaluation reads | selection, column, name |
| `export.nth` | aggregate expression | expression, integer position, default | AnnplyrExpr | caller-defined | requested position | immutable expression | adapter chosen by caller | caller projected budget | position, expression |
| `export.ntile` | row expression | expression and positive bucket count | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | bucket count, expression |
| `export.num_range` | schema-dependent | prefix, integer range, width | tidy selector | caller-defined | requested numeric order | immutable helper | schema-only | no direct read | unknown column, selection |
| `export.obs_names` | obs rows | virtual axis names | AnnplyrExpr with virtual exact dependency | caller-defined | obs positional order | immutable expression | no matrix read | caller budget excludes metadata | incompatible axis |
| `export.percent_rank` | row expression | string or expression | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, size |
| `export.pick` | schema-dependent | tidy selector | expression-selection helper | caller-defined | selector order | immutable helper | schema-only until caller evaluates | caller projected budget | selection |
| `export.recode` | row expression | expression, mapping, optional default | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, dtype |
| `export.rename_obs_names` | obs names | AnnData and mapper | AnnData | ungrouped | stable | same-shape | backed write rejected; sparse unchanged | no matrix reads | duplicate name, backed write |
| `export.rename_var_names` | var names | AnnData and mapper | AnnData | ungrouped | stable | same-shape | backed write rejected; sparse unchanged | no matrix reads | duplicate name, backed write |
| `export.replace_na` | row expression | expression and value | AnnplyrExpr | caller-defined | row order | immutable expression | adapter chosen by caller | caller projected budget | expression, dtype |
| `export.replace_name_suffix` | obs or var names | AnnData, old/new suffix, axis | AnnData | ungrouped | stable | same-shape | backed write rejected; sparse unchanged | no matrix reads | axis, duplicate name, backed write |
| `export.row_number` | rows | virtual row positions | AnnplyrExpr with no source columns | caller-defined | zero-based evaluation order represented as one-based values | immutable expression | no matrix read | caller budget excludes metadata | expression context |
| `export.sample_meta` | obs | AnnData and sample column/selectors | pandas DataFrame | grouped by sample key | first-seen sample order | independent table | metadata only; backed safe | no matrix reads | selection, unknown column |
| `export.sd` | aggregate expression | string or expression | AnnplyrExpr | caller-defined | scalar | immutable expression | adapter chosen by caller | caller projected budget | expression |
| `export.starts_with` | schema-dependent | prefix | tidy selector | caller-defined | schema order | immutable helper | schema-only | no direct read | selection |
| `export.store_palette` | obs/uns | AnnData, categorical column, palette/key | AnnData | ungrouped | category order | same-shape | backed write rejected | no matrix reads | column, palette size, backed write |
| `export.sum` | aggregate expression | string or expression | AnnplyrExpr | caller-defined | scalar | immutable expression | adapter chosen by caller | caller projected budget | expression |
| `export.tidy` | AnnData | AnnData object | TidyAnnData | caller-defined | caller-defined | upgrades the object in place, keeping identity | no data read | no direct read | none |
| `export.unnest` | obs rows | table with a column of AnnData subsets | AnnData | none | table row order, then subset order | independent concatenation | `anndata.concat` with outer join | no matrix evaluation reads | selection |
| `export.var_names` | var rows | virtual axis names | AnnplyrExpr with virtual exact dependency | caller-defined | var positional order | immutable expression | no matrix read | caller budget excludes metadata | incompatible axis |
| `export.where` | schema-dependent | dtype/schema predicate | tidy selector | caller-defined | schema order | immutable helper | zero-length typed Series; no values read | no direct read | value-dependent predicate, selection |
| `export.add_count` | selected metadata axis | metadata grouping selectors and optional weight; adata passed as the first argument | AnnData | ungrouped | stable unless sort requests axis reorder | same-shape | backed write rejected; sparse format preserved | no matrix reads | selection, column, axis, size, backed write |
| `export.aggregate_cells` | obs rows | AnnData, obs key selectors, required source matrix, fun; adata passed as the first argument | AnnData | one row per observed key combination | first-seen key order | new AnnData with keys and `n_cells` in obs | sparse one-hot product; reads the whole source | no projected budget | selection, source, column, axis |
| `export.anti_join` | selected axis | pandas DataFrame or mapping metadata; adata passed as the first argument | AnnData | ungrouped | surviving left positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | column, axis, NA mode |
| `export.arrange` | obs and var | obs, var, X/layer, raw, obsm, varm expressions; adata passed as the first argument | AnnData | ungrouped | stable explicit sort, positional ties | shape-copy | preserve sparse format; backed input returns the view | projected budget | selection, source, column, axis, budget |
| `export.count` | selected metadata axis | metadata grouping selectors and optional weight; adata passed as the first argument | pandas DataFrame | ephemeral by selector | first-seen or descending count when sort | result table independent | metadata only; backed safe | no matrix reads | selection, column, axis |
| `export.distinct` | selected axis | obs, var, X expressions; adata passed as the first argument | AnnData | ungrouped | positional first occurrence | shape-copy | preserve sparse format; backed input returns the view | projected budget | selection, source, column, axis, budget |
| `export.filter` | obs and var | obs, var, X/layer, raw, names, obsm, varm expressions; adata passed as the first argument | AnnData | ungrouped | positional survivor order | shape-copy | preserve sparse format; backed input returns the view | projected budget | selection, source, column, axis, budget, size |
| `export.full_join` | selected axis | pandas DataFrame or mapping metadata; adata passed as the first argument | AnnData | ungrouped | left positional order | shape-copy | preserve sparse format; reject right-only records | no matrix reads | join relationship, right-only, column, name repair, axis |
| `export.group_by` | exactly one of obs or var, given as `axis=` | existing metadata column names or tidy schema selectors, one per positional key; adata passed as the first argument | GroupedAnnData or exact AnnData when no axis supplied | persistent GroupSpec | group plan first-seen order | wrapper shares exact underlying object | no matrix, virtual, computed, or empty grouping | no matrix reads | selection, column, incompatible axis |
| `export.inner_join` | selected axis | pandas DataFrame or mapping metadata; adata passed as the first argument | AnnData | ungrouped | surviving left positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | join relationship, column, name repair, axis |
| `export.left_join` | selected axis | pandas DataFrame or mapping metadata; adata passed as the first argument | AnnData | ungrouped | left positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | join relationship, column, name repair, axis |
| `export.mutate` | obs and var metadata | obs, var, X/layer, raw, obsm, varm expressions and pandas aligned sources; adata passed as the first argument | AnnData | ungrouped | axis order stable; assignments strictly sequential | same-shape | preserve sparse format; backed write rejected | projected budget | selection, source, column, size, budget, backed write |
| `export.pipe` | generic | callable or callable/keyword tuple; adata passed as the first argument | callable-defined | passes exact AnnData | callable-defined | callable-defined | callable-defined | callable-defined | duplicate pipe keyword, invalid callable |
| `export.pivot_longer` | obs rows by features | obs, X/layer, raw; adata passed as the first argument | pandas DataFrame | ungrouped | obs-major then selected feature order | independent long table | only selected dense/sparse/backed values read | projected budget | selection, source, column, name repair, budget |
| `export.pull` | source-dependent | exactly one column from obs, var, X/layer, raw, one obsm/varm/obsp/varp key, or one uns key; adata passed as the first argument | pandas Series | ungrouped | source positional order | independent pandas result | dense, sparse, DataFrame, and backed projections | projected budget | selection, source, column, budget |
| `export.relocate` | obs metadata, var metadata, X columns | schema selectors and before/after anchor; adata passed as the first argument | AnnData | ungrouped | requested columns relocated; X relocation reorders var positionally | same-shape | preserve sparse format; backed write rejected | schema-only; no charge | selection, column, axis, backed write |
| `export.rename` | obs metadata, var metadata, X columns | explicit old-to-new mappings; adata passed as the first argument | AnnData | ungrouped | stable except requested names | same-shape | preserve sparse format; backed write rejected | schema-only; no charge | column, duplicate name, backed write |
| `export.rename_with` | obs metadata, var metadata, X columns | callable plus schema selectors; adata passed as the first argument | AnnData | ungrouped | stable except requested names | same-shape | preserve sparse format; backed write rejected | schema-only; no charge | selection, column, duplicate name, backed write |
| `export.right_join` | selected axis | pandas DataFrame or mapping metadata; adata passed as the first argument | AnnData | ungrouped | right-key order then left order within key | shape-copy | preserve sparse format; reject right-only records | no matrix reads | join relationship, right-only, column, name repair, axis |
| `export.select` | obs metadata, var metadata, X columns | obs, var, X schema selectors; adata passed as the first argument | AnnData | ungrouped | requested column order; positional var order for X | shape-copy | preserve sparse format; backed input returns the view | schema-only; no charge | selection, source, column, duplicate name |
| `export.semi_join` | selected axis | pandas DataFrame or mapping metadata; adata passed as the first argument | AnnData | ungrouped | surviving left positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | column, axis, NA mode |
| `export.slice` | selected axis | integer positions and slices; adata passed as the first argument | AnnData | ungrouped | requested positional order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | selection, axis |
| `export.slice_head` | selected axis | nonnegative n; adata passed as the first argument | AnnData | ungrouped | first n positions | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | selection, axis |
| `export.slice_max` | selected axis | metadata expression and n; adata passed as the first argument | AnnData | ungrouped | stable descending expression order | shape-copy | preserve sparse format; backed input returns the view | no matrix source reads | selection, column, axis |
| `export.slice_min` | selected axis | metadata expression and n; adata passed as the first argument | AnnData | ungrouped | stable ascending expression order | shape-copy | preserve sparse format; backed input returns the view | no matrix source reads | selection, column, axis |
| `export.slice_sample` | selected axis | n or prop, replace, random state; adata passed as the first argument | AnnData | ungrouped | RNG selection order | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | selection, axis |
| `export.slice_tail` | selected axis | nonnegative n; adata passed as the first argument | AnnData | ungrouped | last n positions | shape-copy | preserve sparse format; backed input returns the view | no matrix reads | selection, axis |
| `export.summarise` | obs and var evaluation axes | same as summarize; adata passed as the first argument | pandas DataFrame | optional ephemeral by selector | first-seen key order | alias of summarize | sparse/backed projected reads supported | projected budget | same as summarize |
| `export.summarize` | obs and var evaluation axes | obs, var, X/layer, raw, obsm, varm expressions; adata passed as the first argument | pandas DataFrame | optional ephemeral by selector | first-seen key order | result table independent | sparse/backed projected reads supported | projected budget | selection, source, column, size, budget |
| `export.to_wide` | obs rows | obs, X/layer, raw, obsm, obsp; adata passed as the first argument | pandas DataFrame | ungrouped | obs positional order and requested columns | independent pandas result | wide sparse uses pandas sparse columns; backed projection | projected budget | selection, source, column, name repair, budget |
| `export.tidyr` | module | the 16 pandas-only helpers: chop, drop_na, extract, fill, hoist, nest, pack, pivot_wider, separate, separate_rows, unchop, unite, unnest, unnest_longer, unnest_wider, unpack | module | none | per helper | independent table | pandas only | no matrix reads | column, name repair |
| `export.transmute` | obs and var metadata | obs, var, X/layer, raw, obsm, varm expressions and pandas aligned sources; adata passed as the first argument | AnnData | ungrouped | requested assignment order; sequential dependencies | independent | preserve sparse format; backed input yields in-memory result | projected budget | selection, source, column, size, budget |

## Join details

All joins validate `relationship`, `multiple`, `suffixes`, and `na_matches`
before subsetting. `unmatched="error"` rejects unmatched left records;
`unmatched="drop"` adds no restriction beyond the join kind. Right and full
joins reject every right-only record regardless of `unmatched`, because aligned
AnnData storage cannot be synthesized. Left/full preserve left order;
inner/semi/anti preserve surviving left order; right follows right-key order and
original left order within each key.

## Expression and evaluation details

Assignments are strictly sequential: later assignments see earlier new or
overwritten columns. Only a contiguous set with known independent dependencies,
compatible known widths, and compatible cardinality may share evaluation.
Opaque raw Narwhals expressions evaluate alone, project and charge the full
source, and are never batched. Scalar recycling and pandas extension dtypes are
preserved.

All matrix-reading calls plan every source request and cumulative charge before
the first adapter read. Pandas DataFrames, dense arrays/views, CSR/CSC matrices
and arrays, backed dense datasets, `CSRDataset`, and `CSCDataset` are supported.
Backed fancy indexing may sort positions internally only if requested order is
restored. Sparse wide exports retain pandas sparse columns where supported.
