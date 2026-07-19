from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, cast

import pandas as pd
from anndata import AnnData

from annplyr._errors import (
    AnnplyrError,
    DuplicateNameError,
    IncompatibleAxisError,
    SelectionError,
    SizeMismatchError,
    UnknownColumnError,
    UnknownSourceError,
)
from annplyr._frames import evaluate_select, obs_frame
from annplyr._verbs import _axis, _ensure_not_backed, left_join_adata, summarize_adata


@dataclass(frozen=True)
class FeaturePresence:
    """Feature presence diagnostics returned by feature presence checks."""

    found_features: list[str]
    missing_features: list[str]
    wrong_case_found_features: list[str]
    _features: tuple[str, ...] = field(default_factory=tuple, repr=False)
    _wrong_case_by_feature: Mapping[str, str] = field(default_factory=dict, repr=False)

    @property
    def all_found(self) -> bool:
        """Return whether every requested feature was found exactly."""
        return not self.missing_features

    def to_frame(self) -> pd.DataFrame:
        """Return feature diagnostics as a pandas data frame."""
        features = self._features or (*self.found_features, *self.missing_features)
        found = set(self.found_features)
        rows = [
            {
                "feature": feature,
                "found": feature in found,
                "wrong_case_found": self._wrong_case_by_feature.get(feature, pd.NA),
            }
            for feature in features
        ]
        return pd.DataFrame(rows, columns=["feature", "found", "wrong_case_found"])


def sample_meta(
    adata: AnnData,
    sample: str,
    *,
    include: Any = None,
    exclude: Any = None,
    on_conflict: str = "error",
) -> pd.DataFrame:
    """Return one row of sample-level observation metadata per sample."""
    _check_obs_column(adata, sample)
    if on_conflict not in {"error", "first"}:
        msg = "on_conflict must be 'error' or 'first'"
        raise AnnplyrError(msg)

    selected = _selected_obs_columns(adata, include)
    excluded = set(_selected_obs_columns(adata, exclude)) if exclude is not None else set()
    columns = [sample, *[column for column in selected if column != sample and column not in excluded]]

    obs_table = _obs_table(adata)
    work = pd.DataFrame(obs_table.loc[:, columns]).copy()
    grouped = work.groupby(sample, sort=False, dropna=False)
    if on_conflict == "error":
        conflicts = [
            column
            for column in columns
            if column != sample and any(group[column].nunique(dropna=False) > 1 for _, group in grouped)
        ]
        if conflicts:
            msg = f"Column(s) are not sample-constant for {sample!r}: {', '.join(conflicts)}"
            raise AnnplyrError(msg)

    rows = work.drop_duplicates(subset=[sample], keep="first").reset_index(drop=True)
    return rows.loc[:, columns]


def sample_summary(
    adata: AnnData,
    sample: str,
    *,
    obs: Mapping[str, Any] | None = None,
    x: Mapping[str, Any] | None = None,
    layer: str | None = None,
) -> pd.DataFrame:
    """Summarize observation metadata and selected expression by sample."""
    _check_obs_column(adata, sample)
    return summarize_adata(adata, obs=obs, x=x, by=sample, layer=layer)


def add_sample_meta(
    adata: AnnData,
    meta: pd.DataFrame | Mapping[str, Any],
    *,
    sample: str,
    by: str | None = None,
    copy: bool = True,
    suffixes: tuple[str, str] = ("", "_meta"),
) -> AnnData:
    """Add sample-level metadata to ``adata.obs`` without changing observation order."""
    _check_obs_column(adata, sample)
    by = sample if by is None else by
    metadata = _coerce_sample_metadata(meta)
    if by not in metadata.columns:
        msg = f"Unknown metadata join column: {by}"
        raise UnknownColumnError(msg)
    if by != sample:
        if sample in metadata.columns:
            msg = f"Metadata already contains target sample column {sample!r}"
            raise DuplicateNameError(msg)
        metadata = metadata.rename(columns={by: sample})
    return left_join_adata(adata, metadata, by=sample, axis="obs", suffixes=suffixes, copy=copy)


def feature_present(
    adata: AnnData,
    features: str | Sequence[str],
    *,
    feature_column: str | None = None,
    case_check: bool = True,
) -> FeaturePresence:
    """Check requested features against ``var_names`` or a feature metadata column."""
    requested = [features] if isinstance(features, str) else [str(feature) for feature in features]
    if feature_column is None:
        available = [str(feature) for feature in adata.var_names]
    else:
        if feature_column not in adata.var.columns:
            msg = f"Unknown feature column: {feature_column}"
            raise UnknownColumnError(msg)
        available = [str(feature) for feature in adata.var[feature_column]]

    available_set = set(available)
    lower_to_feature: dict[str, str] = {}
    for feature in available:
        lower_to_feature.setdefault(feature.lower(), feature)

    found: list[str] = []
    missing: list[str] = []
    wrong_case_by_feature: dict[str, str] = {}
    for feature in requested:
        if feature in available_set:
            found.append(feature)
            continue
        missing.append(feature)
        if case_check:
            wrong_case = lower_to_feature.get(feature.lower())
            if wrong_case is not None:
                wrong_case_by_feature[feature] = wrong_case

    return FeaturePresence(
        found_features=found,
        missing_features=missing,
        wrong_case_found_features=list(wrong_case_by_feature.values()),
        _features=tuple(requested),
        _wrong_case_by_feature=wrong_case_by_feature,
    )


def rename_obs_names(adata: AnnData, mapper: Callable[[str], str] | Mapping[str, str], *, copy: bool = True) -> AnnData:
    """Rename observation names while preserving AnnData alignment.

    Mapping inputs follow the package rename convention: ``{new_name: old_name}``.
    """
    return _rename_axis_names(adata, mapper, axis="obs", copy=copy)


def rename_var_names(adata: AnnData, mapper: Callable[[str], str] | Mapping[str, str], *, copy: bool = True) -> AnnData:
    """Rename variable names while preserving AnnData alignment.

    Mapping inputs follow the package rename convention: ``{new_name: old_name}``.
    """
    return _rename_axis_names(adata, mapper, axis="var", copy=copy)


def replace_name_suffix(
    adata: AnnData,
    current_suffix: str,
    new_suffix: str = "",
    *,
    axis: str = "obs",
    copy: bool = True,
) -> AnnData:
    """Replace a suffix on observation or variable names."""
    if current_suffix:
        mapper = lambda name: (
            f"{name.removesuffix(current_suffix)}{new_suffix}" if name.endswith(current_suffix) else name
        )
    else:
        mapper = lambda name: f"{name}{new_suffix}"
    return _rename_axis_names(adata, mapper, axis=axis, copy=copy)


def add_name_prefix(adata: AnnData, prefix: str, *, axis: str = "obs", sep: str = "_", copy: bool = True) -> AnnData:
    """Add a prefix to observation or variable names."""
    return _rename_axis_names(adata, lambda name: f"{prefix}{sep}{name}", axis=axis, copy=copy)


def name_duplicates(adata: AnnData, *, axis: str = "obs") -> pd.DataFrame:
    """Report duplicated observation or variable names."""
    names = list(_axis_names(adata, axis))
    positions_by_name: dict[str, list[int]] = {}
    for position, name in enumerate(names):
        positions_by_name.setdefault(str(name), []).append(position)
    rows = [
        {"name": name, "count": len(positions), "positions": positions}
        for name, positions in positions_by_name.items()
        if len(positions) > 1
    ]
    return pd.DataFrame(rows, columns=["name", "count", "positions"])


def store_palette(
    adata: AnnData,
    obs: str,
    palette: Mapping[Any, str] | Sequence[str],
    *,
    key: str | None = None,
    copy: bool = True,
) -> AnnData:
    """Store a Scanpy-compatible observation palette in ``adata.uns``."""
    _ensure_not_backed(adata, "store_palette")
    levels = _obs_levels(adata, obs)
    colors = _palette_colors(levels, palette)
    out = adata.copy() if copy else adata
    out.uns[key or f"{obs}_colors"] = colors
    return out


def get_palette(adata: AnnData, obs: str, *, key: str | None = None) -> dict[Any, str]:
    """Return a Scanpy-compatible observation palette as a level-to-color mapping."""
    levels = _obs_levels(adata, obs)
    palette_key = key or f"{obs}_colors"
    if palette_key not in adata.uns:
        msg = f"Unknown palette key: {palette_key}"
        raise UnknownSourceError(msg)
    colors = list(adata.uns[palette_key])
    if len(colors) != len(levels):
        msg = f"Palette length {len(colors)} does not match {len(levels)} level(s) for {obs!r}"
        raise SizeMismatchError(msg)
    return dict(zip(levels, colors, strict=True))


def _selected_obs_columns(adata: AnnData, selector: Any) -> list[str]:
    if selector is None:
        return [str(column) for column in _obs_table(adata).columns]
    selected = evaluate_select(obs_frame(adata), selector)
    obs_columns = _obs_table(adata).columns
    selected_columns = [str(column) for column in selected.columns]
    computed = [column for column in selected_columns if column not in obs_columns]
    if computed:
        msg = f"sample metadata selectors must resolve to existing obs metadata columns: {', '.join(computed)}"
        raise SelectionError(msg)
    return selected_columns


def _coerce_sample_metadata(meta: pd.DataFrame | Mapping[str, Any]) -> pd.DataFrame:
    try:
        return pd.DataFrame(meta).copy()
    except ValueError as exc:
        msg = "sample metadata must be a pandas DataFrame or mapping/list of column vectors or records"
        raise AnnplyrError(msg) from exc


def _check_obs_column(adata: AnnData, column: str) -> None:
    if column not in _obs_table(adata).columns:
        msg = f"Unknown obs column: {column}"
        raise UnknownColumnError(msg)


def _rename_axis_names(
    adata: AnnData,
    mapper: Callable[[str], str] | Mapping[str, str],
    *,
    axis: str,
    copy: bool,
) -> AnnData:
    _ensure_not_backed(adata, f"rename_{axis}_names")
    axis = _axis(axis)
    names = [str(name) for name in _axis_names(adata, axis)]
    new_names = _map_names(names, mapper, axis=axis)
    _check_unique_names(new_names, axis=axis)
    out = adata.copy() if copy else adata
    if axis == "obs":
        out.obs_names = new_names
    else:
        out.var_names = new_names
    return out


def _map_names(names: Sequence[str], mapper: Callable[[str], str] | Mapping[str, str], *, axis: str) -> list[str]:
    if callable(mapper):
        mapped = [str(mapper(name)) for name in names]
    elif isinstance(mapper, Mapping):
        old_names = list(mapper.values())
        duplicated_sources = pd.Index(old_names)[pd.Index(old_names).duplicated()].unique().tolist()
        if duplicated_sources:
            msg = f"Duplicate {axis} source name(s): {', '.join(str(name) for name in duplicated_sources)}"
            raise DuplicateNameError(msg)
        missing = [old for old in old_names if old not in names]
        if missing:
            msg = f"Unknown {axis} name(s): {', '.join(str(name) for name in missing)}"
            raise UnknownColumnError(msg)
        old_to_new = {str(old): str(new) for new, old in mapper.items()}
        mapped = [old_to_new.get(name, name) for name in names]
    else:
        msg = "name mapper must be a callable or mapping of new names to old names"
        raise AnnplyrError(msg)
    if len(mapped) != len(names):
        msg = f"Renaming {axis} names must preserve length"
        raise SizeMismatchError(msg)
    return mapped


def _check_unique_names(names: Sequence[str], *, axis: str) -> None:
    duplicated = pd.Index(names)[pd.Index(names).duplicated()].unique().tolist()
    if duplicated:
        msg = f"Duplicate {axis} name(s) after operation: {', '.join(str(name) for name in duplicated)}"
        raise DuplicateNameError(msg)


def _axis_names(adata: AnnData, axis: str) -> pd.Index:
    axis = _axis(axis)
    if axis == "obs":
        return adata.obs_names
    if axis == "var":
        return adata.var_names
    raise IncompatibleAxisError("axis must be 'obs' or 'var'")


def _obs_levels(adata: AnnData, obs: str) -> list[Any]:
    _check_obs_column(adata, obs)
    values = _obs_table(adata)[obs]
    if isinstance(values.dtype, pd.CategoricalDtype):
        return list(values.cat.categories)
    return list(pd.unique(values))


def _obs_table(adata: AnnData) -> pd.DataFrame:
    return cast(pd.DataFrame, adata.obs)


def _palette_colors(levels: Sequence[Any], palette: Mapping[Any, str] | Sequence[str]) -> list[str]:
    if isinstance(palette, Mapping):
        missing = [level for level in levels if level not in palette]
        if missing:
            msg = f"Palette is missing level(s): {', '.join(str(level) for level in missing)}"
            raise UnknownColumnError(msg)
        return [str(palette[level]) for level in levels]
    colors = [str(color) for color in palette]
    if len(colors) != len(levels):
        msg = f"Palette length {len(colors)} does not match {len(levels)} level(s)"
        raise SizeMismatchError(msg)
    return colors
