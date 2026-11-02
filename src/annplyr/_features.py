"""``join_features``: bring matrix features next to the cell metadata.

This is tidySingleCellExperiment's ``join_features(features, shape=...)``.
``shape="long"`` returns a DataFrame with one row per cell and feature, with
every obs column plus ``.cell``, ``.feature`` and one ``.abundance_<source>``
column per source. ``shape="wide"`` returns the object with each feature added
to ``obs`` as a column named after the feature.
"""

from __future__ import annotations

from typing import Any, Literal, cast

import pandas as pd
from anndata import AnnData

from annplyr._errors import DuplicateNameError, IncompatibleAxisError, SelectionError
from annplyr._expr import AnnplyrExpr, col, everything
from annplyr._resolve import resolve_source
from annplyr._verbs import mutate_adata, to_tidy_adata

_MATRIX_KINDS = {"x", "raw", "layer"}


def _feature_text(feature: Any) -> str:
    if isinstance(feature, str):
        return feature
    name = getattr(feature, "_source_column", None) if isinstance(feature, AnnplyrExpr) else None
    if name is None:
        raise SelectionError("join_features takes feature names such as raw.MS4A1 or 'raw:MS4A1', not expressions")
    return name


def _by_source(adata: AnnData, features: tuple[Any, ...]) -> dict[tuple[str, str | None], list[str]]:
    """Group the requested features by the matrix they come from, in order."""
    if not features:
        raise SelectionError("join_features needs at least one feature, for example raw.MS4A1")
    groups: dict[tuple[str, str | None], list[str]] = {}
    for feature in features:
        text = _feature_text(feature)
        resolved = resolve_source(adata, text)
        if resolved.kind not in _MATRIX_KINDS:
            bare = resolved.bare_names[0]
            msg = f"join_features reads X, raw or a layer; {text!r} is {resolved.context}. Name the matrix, e.g. raw.{bare}"
            raise IncompatibleAxisError(msg)
        groups.setdefault((resolved.kind, resolved.key), []).extend(resolved.bare_names)
    return groups


def _label(kind: str, key: str | None) -> str:
    return key if kind == "layer" and key is not None else kind


def join_features_adata(
    adata: AnnData,
    *features: Any,
    shape: Literal["long", "wide"] = "long",
    max_matrix_values: int | None = None,
) -> Any:
    groups = _by_source(adata, features)
    if shape == "wide":
        names = [_feature_text(feature) for feature in features]
        bare = [name.partition(":")[2] if ":" in name else name for name in names]
        clashes = sorted({name for name in bare if bare.count(name) > 1} | (set(bare) & set(adata.obs.columns)))
        if clashes:
            msg = f"join_features(shape='wide') would write {clashes} twice; name them with mutate(new=raw.GENE)"
            raise DuplicateNameError(msg)
        return mutate_adata(
            adata, {b: col(n) for b, n in zip(bare, names, strict=True)}, max_matrix_values=max_matrix_values
        )
    if shape != "long":
        raise SelectionError(f"join_features shape= must be 'long' or 'wide', not {shape!r}")

    # Assemble by position, never by name: duplicate cell names stay distinct.
    if len({tuple(names) for names in groups.values()}) > 1:
        raise SelectionError("join_features(shape='long') needs the same features from every source")
    long: pd.DataFrame | None = None
    for (kind, key), names in groups.items():
        value = f".abundance_{_label(kind, key)}"
        part = to_tidy_adata(
            adata,
            obs=everything() if long is None else None,
            raw=names if kind == "raw" else None,
            x=names if kind != "raw" else None,
            layer=key if kind == "layer" else None,
            obs_name=".cell",
            feature=".feature",
            value=value,
            max_matrix_values=max_matrix_values,
            context="join_features",
        )
        if long is None:
            long = part
            continue
        long[value] = part[value].to_numpy()
    return cast(pd.DataFrame, long)
