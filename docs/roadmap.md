# Project direction

`annplyr` aims to be the small, dependable wrangling layer between AnnData and
the wider Python dataframe ecosystem. Development is guided by real single-cell
workflows rather than by matching every dataframe method or adding biological
analysis wrappers.

## Design principles

- **Preserve AnnData alignment.** Every AnnData-returning operation must keep
  `obs`, `var`, `X`, layers, `raw`, `obsm`, `varm`, `obsp`, and `varp`
  coordinated.
- **Keep ownership explicit.** Independent results, optional non-copy results,
  in-place metadata updates, and backed behavior must remain predictable.
- **Project before converting.** Sparse and backed inputs should read selected
  rows and columns before pandas conversion, with no silent whole-matrix
  densification.
- **Prefer composable verbs.** A small set of consistent operations is more
  valuable than aliases, special cases, and overlapping helpers.
- **Fail on unsafe tidy semantics.** A join or reshape that would manufacture
  cells, features, or unaligned records should raise a clear package error.
- **Measure performance claims.** Optimizations must preserve exact results and
  be supported by reproducible timing and memory evidence.

## Current scope

The supported surface includes:

- filtering, selecting, ordering, slicing, mutation, transmutation, summaries,
  and counts for observation and feature axes;
- persistent grouping with deterministic categorical and missing-value
  semantics;
- metadata joins with cardinality and unmatched-row checks;
- expressions and tidyselect-style selectors across metadata, selected `X` or
  layer values, `raw`, embeddings, and loadings;
- wide, long, nested, and controlled container extraction to pandas;
- dense, CSR, CSC, view, and backed source adapters with positional projection;
- cumulative materialization budgets and bounded canonical reductions;
- consumer typing for accessor, grouped, expression, selector, source, and join
  workflows;
- narrow utilities for sample metadata, feature presence, safe names, and
  Scanpy-compatible palettes.

## Near-term priorities

Future work should begin with a demonstrated workflow and remain independently
testable. Likely priorities are:

1. broaden downstream integration fixtures when a real consumer exposes a new
   alignment, dtype, serialization, or backed-storage contract;
2. improve diagnostics where a valid operation fails without enough source,
   axis, or projection context;
3. deepen typing only where downstream projects encounter a concrete gap;
4. optimize measured bottlenecks without adding public tuning controls or
   weakening exact eager equivalence;
5. keep examples representative of current Scanpy and AnnData practice.

Repository issues define accepted work; this page does not promise an API or
timeline for speculative ideas.

## Deliberate non-goals

- biological QC doctrine, marker discovery, cluster annotation, or curated
  species/pathway gene sets;
- plotting wrappers tied to one visualization library;
- direct replacements for Scanpy, scvi-tools, CellRank, or AnnData storage;
- lossy joins that duplicate cells or features;
- whole-matrix long exports by default;
- public chunk-size knobs or backend-specific execution flags;
- full lazy execution across arbitrary dataframe and array backends.

These boundaries keep annplyr focused: make aligned AnnData wrangling concise,
safe, and easy to compose with the scientific tools that come before and after
it.
