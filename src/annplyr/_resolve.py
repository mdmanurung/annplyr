"""Resolve ``source:name`` column references into concrete AnnData sources.

Expressions name their own source, so a verb takes one flat payload instead of
one keyword per storage container. The rule is deliberately narrow:

- Metadata is the default namespace. ``ap.col("n_genes")`` reads ``obs``.
- A matrix always needs a qualifier, because there is no safe default:
  ``x:``, ``raw:``, a layer name such as ``counts:``, or an ``obsm``/``varm``
  key such as ``X_pca:``.
- ``var:``/``varm:`` reference the column axis, so those verbs subset features.

Resolution never rewrites an expression. The projection is keyed by the
qualified reference, so ``ap.col("raw:LYZ")`` is evaluated verbatim against a
frame whose column is literally named ``raw:LYZ``. That keeps the Narwhals
expression identical to what the caller wrote and avoids tree rewriting.

Qualifier parsing is split deliberately. Reserved prefixes are recognized here
with no AnnData access, so ``ENSG00000123456:1``-style feature names are never
mangled. Layer and ``obsm``/``varm`` prefixes need the object, so
:func:`resolve_source` looks them up only after the reserved check fails.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import pandas as pd
from anndata import AnnData

from annplyr._errors import IncompatibleAxisError, UnknownSourceError
from annplyr._expr import col, collect_column_references
from annplyr._frames import OBS_NAMES, VAR_NAMES, obs_frame, var_frame
from annplyr._sources import RequestPlanner

#: Prefixes that mean a specific container without needing the AnnData object.
#: The virtual axis-name columns are deliberately absent: they are already
#: ordinary columns of the obs/var frames, so they need no special handling.
RESERVED_QUALIFIERS = frozenset({"obs", "var", "x", "raw", "obsm", "varm", "obsp", "varp", "uns"})

#: The virtual axis-name columns are ordinary columns of the obs/var frames,
#: but they must still resolve to the right axis, so they are named here rather
#: than falling through to the bare-name ``obs`` default.
_VIRTUAL_AXIS = {OBS_NAMES: "obs", VAR_NAMES: "var"}

#: Which AnnData axis each source subsetting touches.
_AXIS_BY_KIND = {
    "obs": "obs",
    "var": "var",
    "x": "obs",
    "raw": "obs",
    "layer": "obs",
    "obsm": "obs",
    "varm": "var",
    "obsp": "obs",
    "varp": "var",
    "uns": "obs",
}


@dataclass(frozen=True)
class ResolvedSource:
    """One concrete container an expression can be evaluated against."""

    kind: str
    key: str | None
    axis: str
    references: tuple[tuple[str, str], ...]

    @property
    def context(self) -> str:
        return f"{self.kind}:{self.key}" if self.key else self.kind

    @property
    def bare_names(self) -> tuple[str, ...]:
        return tuple(bare for _, bare in self.references)

    @property
    def is_matrix(self) -> bool:
        return self.kind not in {"obs", "var", "names"}


def split_qualifier(text: str) -> tuple[str | None, str]:
    """Split a reserved ``source:name`` prefix, leaving the name otherwise."""
    head, separator, tail = text.partition(":")
    if separator and tail and head in RESERVED_QUALIFIERS:
        return head, tail
    return None, text


def resolve_source(adata: AnnData, text: str) -> ResolvedSource:
    """Resolve one column reference against ``adata``.

    Raises :class:`UnknownSourceError` for an unknown layer, key, or raw
    matrix, and for a colon-prefixed text that matches no container and no
    literal column name, which is almost always a mistyped qualifier.
    """
    qualifier, name = split_qualifier(text)
    if qualifier is None:
        qualifier = _dynamic_qualifier(adata, text)
        if qualifier is not None:
            name = text.partition(":")[2]
    kind, key = _classify(adata, qualifier, text)
    if kind is None:
        raise UnknownSourceError(_unknown_message(adata, text))
    if kind == "obs" and qualifier is None and _looks_mistyped(adata, text):
        raise UnknownSourceError(_unknown_message(adata, text))
    return ResolvedSource(kind, key, _AXIS_BY_KIND[kind], ((text, name),))


def _dynamic_qualifier(adata: AnnData, text: str) -> str | None:
    head, separator, tail = text.partition(":")
    if not separator or not tail:
        return None
    if head in adata.layers or head in adata.obsm or head in adata.varm:
        return head
    return None


def _classify(adata: AnnData, qualifier: str | None, text: str) -> tuple[str | None, str | None]:
    if qualifier is None:
        return _VIRTUAL_AXIS.get(text, "obs"), None
    if qualifier == "obs":
        return "obs", None
    if qualifier == "var":
        return "var", None
    if qualifier == "x":
        return "x", None
    if qualifier == "raw":
        return "raw", None
    if qualifier in {"obsm", "varm", "obsp", "varp", "uns"}:
        return qualifier, text.partition(":")[0]
    if qualifier in adata.layers:
        return "layer", qualifier
    if qualifier in adata.obsm:
        return "obsm", qualifier
    if qualifier in adata.varm:
        return "varm", qualifier
    return None, None


def _looks_mistyped(adata: AnnData, text: str) -> bool:
    """Report whether a colon-prefixed text is probably a mistyped qualifier.

    A feature or cell name may legitimately contain a colon (``ENSG001:1``), so
    a text whose prefix matches no container is only reported when nothing in
    the object actually carries that literal name.
    """
    head, separator, tail = text.partition(":")
    if not separator or not tail:
        return False
    if head in adata.layers or head in adata.obsm or head in adata.varm:
        return False
    if head in RESERVED_QUALIFIERS:
        return False
    known = text in adata.obs.columns or text in adata.var.columns or text in adata.obs_names or text in adata.var_names
    if adata.raw is not None and text in adata.raw.var_names:
        known = True
    return not known


def _unknown_message(adata: AnnData, text: str) -> str:
    known: list[str] = sorted(RESERVED_QUALIFIERS)
    known += [f"{name}:" for name in adata.layers.keys() if name is not None]  # anndata 0.13 lists X as None
    known += [f"{name}:" for name in adata.obsm]
    known += [f"{name}:" for name in adata.varm]
    head = text.partition(":")[0]
    return (
        f"Cannot resolve column reference {text!r}: {head!r} is not a known container. "
        f"Prefix a matrix column with its container, for example 'raw:LYZ', or use a bare "
        f"name for obs metadata. Known prefixes: {', '.join(sorted(known))}"
    )


def collect_references(values: Any) -> set[str] | None:
    """Return every column reference in ``values``, or ``None`` if opaque."""
    return collect_column_references(values)


def group_by_source(adata: AnnData, references: set[str]) -> dict[ResolvedSource, list[str]]:
    """Partition references into one :class:`ResolvedSource` per container.

    The returned keys carry every reference that targets that container, so a
    caller can project exactly those columns.
    """
    grouped: dict[tuple[str, str | None], list[tuple[str, str]]] = {}
    for text in sorted(references):
        source = resolve_source(adata, text)
        grouped.setdefault((source.kind, source.key), []).extend(source.references)
    sources: dict[ResolvedSource, list[str]] = {}
    for (kind, key), pairs in grouped.items():
        source = ResolvedSource(kind, key, _AXIS_BY_KIND[kind], tuple(pairs))
        sources.setdefault(source, []).extend(text for text, _ in pairs)
    return sources


def assert_single_source(adata: AnnData, values: Any) -> None:
    """Reject an expression that spans more than one container.

    Cross-source row-level arithmetic is not supported today either, because
    each source is evaluated against its own frame. Failing here gives a clear
    message instead of a confusing missing-column error later.
    """
    references = collect_references(values)
    if not references:
        return
    sources = group_by_source(adata, references)
    if len(sources) > 1:
        listed = ", ".join(sorted(source.context for source in sources))
        msg = (
            f"one expression cannot mix columns from {listed}. Split it into one "
            f"expression per source, or precompute the shared part with mutate()"
        )
        raise IncompatibleAxisError(msg)


def plan_source(
    planner: RequestPlanner,
    adata: AnnData,
    source: ResolvedSource,
    *,
    row_positions: Any = None,
) -> int:
    """Register one container read and return its token.

    Callers must call ``planner.execute()`` once after planning every source so
    the cumulative budget is validated before the first adapter read.
    """
    from annplyr._sources import source_adapter

    # A layer is the X matrix read through the layer selector, so it is
    # requested as an x adapter with a layer name rather than its own kind.
    adapter_kind = "x" if source.kind == "layer" else source.kind
    layer = source.key if source.kind == "layer" else None
    key = None if source.kind == "layer" else source.key
    adapter = source_adapter(adata, adapter_kind, key=key, layer=layer)
    # An expression that named no column, such as a raw Narwhals expression,
    # gets a None request so the planner charges the whole container rather than
    # projecting zero columns.
    request = col(*source.bare_names) if source.bare_names else None
    return planner.add(
        adapter,
        request=request,
        mode="expression",
        row_positions=row_positions,
        context=source.context,
        source=source,
    )


def resolve_container(adata: AnnData, name: str) -> ResolvedSource:
    """Resolve a bare container name such as ``raw``, ``x``, or ``counts``.

    Used for expressions that name no column, such as a raw Narwhals
    expression, which cannot declare its own container. Every such expression is
    evaluated against the whole container, so the planner charges it in full.
    """
    if name in adata.layers:
        kind, key = "layer", name
    elif name in adata.obsm:
        kind, key = "obsm", name
    elif name in adata.varm:
        kind, key = "varm", name
    elif name in {"obs", "var", "x", "raw"}:
        kind, key = name, None
    else:
        raise UnknownSourceError(_unknown_message(adata, f"{name}:column"))
    return ResolvedSource(kind, key, _AXIS_BY_KIND[kind], ())


def partition_payload(
    adata: AnnData,
    items: Sequence[Any],
    *,
    verb: str,
    default_source: ResolvedSource | None = None,
) -> dict[ResolvedSource, list[Any]]:
    """Split individual payload items by the one container each one references.

    An item that references no column at all — a literal, or an opaque Narwhals
    expression whose references are unknowable — is evaluated against
    ``default_source``, or the ``obs`` frame when none is given. An item spanning
    two containers is rejected, because each container is evaluated against its
    own frame and cannot be combined at row level.
    """
    buckets: dict[ResolvedSource, list[Any]] = {}
    leftovers: list[Any] = []
    for item in items:
        references = collect_column_references(item)
        if not references:
            leftovers.append(item)
            continue
        sources = group_by_source(adata, references)
        if len(sources) > 1:
            listed = ", ".join(sorted(source.context for source in sources))
            msg = (
                f"{verb}: one expression cannot mix columns from {listed}. "
                f"Split it into one expression per container, or precompute the shared "
                f"part with mutate() first"
            )
            raise IncompatibleAxisError(msg)
        source = next(iter(sources))
        buckets.setdefault(source, []).append(item)
    if leftovers:
        fallback = default_source if default_source is not None else ResolvedSource("obs", None, "obs", ())
        buckets.setdefault(fallback, []).extend(leftovers)
    return buckets


def prepare_frames(
    adata: AnnData,
    buckets: Mapping[ResolvedSource, Sequence[Any]],
    planner: RequestPlanner,
    *,
    row_positions: Any = None,
    materialize: bool = True,
) -> dict[ResolvedSource, pd.DataFrame]:
    """Plan, validate, and read every container a payload needs.

    All matrix reads are registered before any adapter read, so the cumulative
    budget is validated first. Returned frames carry the caller's reference text
    as their column names, so expressions evaluate verbatim.

    With ``materialize=False`` only the plan and its budget check run, and matrix
    containers are omitted from the result. That is for callers that read a
    container themselves in bounded chunks, so the projection is read once.
    """
    tokens = {
        source: plan_source(planner, adata, source, row_positions=row_positions)
        for source in buckets
        if source.is_matrix
    }
    if materialize:
        frames = _read_frames(adata, buckets, planner, tokens)
    else:
        # validate the cumulative budget without reading, for callers that read
        # the container themselves in bounded chunks
        planner.validate()
        frames = {}
    return frames


def _read_frames(
    adata: AnnData,
    buckets: Mapping[ResolvedSource, Sequence[Any]],
    planner: RequestPlanner,
    tokens: Mapping[ResolvedSource, int],
) -> dict[ResolvedSource, pd.DataFrame]:
    projected = planner.execute() if tokens else []
    frames: dict[ResolvedSource, pd.DataFrame] = {}
    for source in buckets:
        if source.is_matrix:
            frames[source] = label_frame(projected[tokens[source]], source)
        else:
            frames[source] = metadata_frame(adata, source)
    return frames


def planner_tokens(
    planner: RequestPlanner,
    adata: AnnData,
    buckets: Mapping[ResolvedSource, Sequence[Any]],
    *,
    row_positions: Any = None,
) -> dict[ResolvedSource, int]:
    """Register every matrix read and return its token per source."""
    return {
        source: plan_source(planner, adata, source, row_positions=row_positions)
        for source in buckets
        if source.is_matrix
    }


def label_frame(frame: pd.DataFrame, source: ResolvedSource) -> pd.DataFrame:
    """Rename a projected frame's columns to their qualified reference text.

    ``references`` holds ``(reference_text, bare_name)`` pairs, so the rename
    maps each bare column name onto the text the caller wrote.
    """
    return frame.rename(columns={bare: reference for reference, bare in source.references})


def metadata_frame(adata: AnnData, source: ResolvedSource) -> pd.DataFrame:
    """Read a metadata container, naming columns by qualified reference."""
    frame = obs_frame(adata) if source.kind == "obs" else var_frame(adata)
    return label_frame(frame, source)


def axis_sources(sources: Sequence[ResolvedSource]) -> tuple[list[ResolvedSource], list[ResolvedSource]]:
    """Split sources into row-axis and column-axis halves."""
    rows = [source for source in sources if source.axis == "obs"]
    columns = [source for source in sources if source.axis == "var"]
    return rows, columns


def resolve_grouping_keys(selector: Any, *, axis: str | None = None) -> Any:
    """Validate grouping keys against the axis being reduced, returning bare names.

    A grouping key is metadata on the reduced axis. How much it has to say about
    that depends on whether the verb already knows:

    - ``axis=None`` (``group_by``, which takes an explicit ``axis=``): the axis is
      known, so a qualifier is redundant and is rejected.
    - ``axis="obs"`` or ``"var"`` (``summarize``, ``count``, which learn the axis
      from the aggregates): a bare name is ``obs`` and a ``var:`` name is a
      feature, so the qualifier is required to disambiguate a feature key.

    A key naming the axis other than the one being reduced is rejected rather
    than resolved. Resolving ``by="var:kind"`` against an ``obs`` frame would
    silently group whichever ``obs`` column happens to share the bare name, which
    is the wrong answer rather than an error.
    """
    if isinstance(selector, str):
        qualifier, separator, tail = selector.partition(":")
        if separator and tail and qualifier in RESERVED_QUALIFIERS:
            if axis is None:
                raise IncompatibleAxisError(
                    f"the axis is already given, so {selector!r} is redundant: pass the bare metadata name instead"
                )
            if qualifier not in (axis, "var" if axis == "var" else "obs"):
                raise IncompatibleAxisError(
                    f"{selector!r} names the {qualifier} axis, but this call reduces the {axis} "
                    f"axis: group by a bare name for {axis} metadata, or {qualifier}-qualify the "
                    "aggregates to reduce that axis"
                )
            return tail
        if axis == "var":
            raise IncompatibleAxisError(
                f"a bare name is obs metadata, so grouping the var axis needs a qualifier: use {selector!r}"
            )
        return selector
    if isinstance(selector, Mapping):
        return {key: resolve_grouping_keys(value, axis=axis) for key, value in selector.items()}
    if isinstance(selector, Sequence) and not isinstance(selector, (str, bytes)):
        return [resolve_grouping_keys(item, axis=axis) for item in selector]
    return selector
