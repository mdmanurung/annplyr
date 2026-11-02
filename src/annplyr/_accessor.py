from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any, TypeVar, cast, overload

import pandas as pd
from anndata import AnnData

from annplyr._errors import AnnplyrError, IncompatibleAxisError
from annplyr._expr import AnnplyrExpr, everything
from annplyr._extensions import register_anndata_accessor
from annplyr._grouped import GroupedAnnData
from annplyr._resolve import (
    META_AND_X,
    container_kwargs,
    key_names,
    relocate_anchor,
    rename_kwargs,
)
from annplyr._tidy import TidyAnnData, install_verbs
from annplyr._typing import (
    AnnDataWithAnnplyr,
    Axis,
    Expression,
    JoinBy,
    JoinInput,
    JoinMultiple,
    JoinRelationship,
    JoinUnmatched,
    NaMatches,
    Selector,
    Source,
)
from annplyr._verbs import (
    _validate_filter_scope,
    add_count_adata,
    add_tally_adata,
    anti_join_adata,
    arrange_adata,
    as_frame_adata,
    count_adata,
    distinct_adata,
    filter_adata,
    full_join_adata,
    inner_join_adata,
    left_join_adata,
    merge_named,
    mutate_adata,
    nest_by_adata,
    pivot_longer_adata,
    pull_adata,
    relocate_adata,
    rename_adata,
    rename_with_adata,
    right_join_adata,
    select_adata,
    semi_join_adata,
    slice_adata,
    slice_head_adata,
    slice_max_adata,
    slice_min_adata,
    slice_sample_adata,
    slice_tail_adata,
    summarize_adata,
    tally_adata,
    to_tidy_adata,
    to_wide_adata,
    transmute_adata,
)

T = TypeVar("T")

_TABLE = frozenset({"obs", "x", "raw", "layer", "obsm", "obsp"})

_PULLABLE = _TABLE | {"var", "varm", "varp", "uns"}


def _with_annplyr_accessor(adata: AnnData) -> AnnDataWithAnnplyr:
    """Return a chainable result that still passes ``isinstance(x, AnnData)``.

    Wrapping preserves object identity, so ``inplace=True`` verbs still return
    the exact object they were called on.
    """
    return cast(AnnDataWithAnnplyr, TidyAnnData.wrap(adata))


@register_anndata_accessor("ap")
class AnnplyrAccessor:
    """Dataframe-style AnnData wrangling accessor."""

    def __init__(self, adata: AnnData):
        self._obj = cast(AnnDataWithAnnplyr, adata)

    def filter(
        self,
        *predicates: Expression | Sequence[Expression],
        source: str | None = None,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        """Subset axes with one predicate per argument, combined with AND.

        Each predicate names its own container, so a single call can read
        ``obs`` metadata, ``raw``, a named layer, and an embedding together::

            adata.ap.filter(
                ap.col("n_genes") >= 500,
                ap.col("raw:LYZ") < 1.0,
                ap.col("counts:MS4A1") > 0,
            )

        A bare name is ``obs`` metadata. A matrix always needs a qualifier,
        because there is no safe default container for it: use ``x:`` for ``X``,
        ``raw:`` for ``.raw``, a layer name such as ``counts:``, or an
        ``obsm``/``varm`` key such as ``X_pca:``. ``var:`` and ``varm:``
        reference the feature axis and subset features.

        One expression cannot mix containers, because each is evaluated against
        its own frame. Precompute the shared part with :meth:`mutate` first.

        ``source=`` names the container for expressions that reference no
        column, which is how a raw Narwhals expression is pointed at a matrix.
        Such an expression is charged for the whole container.

        Parameters
        ----------
        *predicates
            See the annotated public signature.
        copy
            See the annotated public signature.
        max_matrix_values
            See the annotated public signature.

        Returns
        -------
        AnnData

        Ownership
        ---------
        Returns an independent result by default; ``copy=False`` may be a view or materialized.

        Raises
        ------
        SelectionError
            If the request violates the corresponding public contract.
        UnknownColumnError
            If the request violates the corresponding public contract.
        UnknownSourceError
            If the request violates the corresponding public contract.
        IncompatibleAxisError
            If the request violates the corresponding public contract.
        AnnplyrError
            If the request violates the corresponding public contract.
        """
        return _with_annplyr_accessor(
            filter_adata(
                self._obj,
                *predicates,
                source=source,
                copy=copy,
                max_matrix_values=max_matrix_values,
            )
        )

    def filter_obs(
        self,
        *predicates: Expression | Sequence[Expression],
        source: str | None = None,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        """Filter cells using only obs metadata, combining predicates with AND.

        Bare names and ``obs:`` references read observation metadata. Sources
        on the feature axis or in matrices are rejected before reads. Ownership,
        null handling, and predicate evaluation follow :meth:`filter`, including
        its no-op behavior when no predicates are supplied.

        Parameters
        ----------
        *predicates
            Boolean expressions, combined with AND.
        source
            Container for an expression that names no column.
        copy
            Return an independent object rather than a view.
        max_matrix_values
            Upper bound on matrix values read; ``None`` means no bound.

        Returns
        -------
        AnnData
            The observations where every predicate is true.

        Ownership
        ---------
        Returns an independent result by default; ``copy=False`` may be a view.

        Raises
        ------
        IncompatibleAxisError
            If a predicate reads a source outside this verb's scope.
        UnknownColumnError
            If a predicate names a column or feature that does not exist.
        """
        _validate_filter_scope(self._obj, predicates, source=source, expression=False)
        return self.filter(*predicates, source=source, copy=copy, max_matrix_values=max_matrix_values)

    def filter_exp(
        self,
        *predicates: Expression | Sequence[Expression],
        source: str | None = None,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        """Filter cells using explicitly qualified gene-expression predicates.

        Use ``ap.col("x:MS4A1")``, ``ap.col("raw:MS4A1")``, or a named layer
        such as ``ap.col("counts:MS4A1")``. Bare names never guess an expression
        source. An opaque Narwhals expression needs ``source="x"``, ``"raw"``,
        or a layer name, and charges the entire source to the matrix budget.
        Metadata, embeddings, and feature-axis sources are rejected before
        reads. Evaluation and ownership follow :meth:`filter`, including its
        no-op behavior when no predicates are supplied.

        Parameters
        ----------
        *predicates
            Boolean expressions, combined with AND.
        source
            Container for an expression that names no column.
        copy
            Return an independent object rather than a view.
        max_matrix_values
            Upper bound on matrix values read; ``None`` means no bound.

        Returns
        -------
        AnnData
            The observations where every predicate is true.

        Ownership
        ---------
        Returns an independent result by default; ``copy=False`` may be a view.

        Raises
        ------
        IncompatibleAxisError
            If a predicate reads a source outside this verb's scope.
        UnknownColumnError
            If a predicate names a column or feature that does not exist.
        """
        _validate_filter_scope(self._obj, predicates, source=source, expression=True)
        return self.filter(*predicates, source=source, copy=copy, max_matrix_values=max_matrix_values)

    def select(self, *names: Selector, copy: bool = True) -> AnnDataWithAnnplyr:
        """Keep the named columns and features; ``-obs.x`` drops one instead.

        Names are strings, container columns (``obs.louvain``, ``var.gene_ids``,
        ``X.MS4A1``) or selectors (``ap.starts_with("n_")``, read as ``obs``).
        Each container is reshaped on its own, so ``select(obs.louvain, X.MS4A1)``
        keeps one obs column and one feature. A whole ``Source`` such as ``X``
        keeps every column of it. With no names nothing changes.
        """
        kwargs = container_kwargs(self._obj, names, verb="select", allowed=META_AND_X, negation=True)
        return _with_annplyr_accessor(select_adata(self._obj, **kwargs, copy=copy))

    def rename(
        self,
        mapping: Mapping[str, Any] | None = None,
        /,
        *,
        inplace: bool = False,
        **new_to_old: Any,
    ) -> AnnDataWithAnnplyr:
        """Rename columns or features, as ``new=old`` where ``old`` is a string or one column.

        ``rename(cell_type=obs.louvain)`` renames an obs column and
        ``rename({"CD4+": X.CD4})`` takes a name that is not an identifier. The
        container comes from the old name, so a bare string is ``obs``. A new name
        that clashes with ``inplace`` goes in the mapping.
        """
        return _with_annplyr_accessor(
            rename_adata(self._obj, **rename_kwargs(self._obj, mapping, new_to_old), inplace=inplace)
        )

    def rename_with(
        self,
        func: Callable[[str], str],
        *names: Selector,
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        """Rename the named columns with ``func``; with no names, every ``obs`` column."""
        kwargs = container_kwargs(self._obj, names, verb="rename_with", allowed=META_AND_X) or {"obs": everything()}
        return _with_annplyr_accessor(rename_with_adata(self._obj, func, **kwargs, inplace=inplace))

    def relocate(
        self,
        *names: Selector,
        before: str | AnnplyrExpr | None = None,
        after: str | AnnplyrExpr | None = None,
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        """Move the named columns to the front, or before or after an anchor column."""
        kwargs = container_kwargs(self._obj, names, verb="relocate", allowed=META_AND_X, negation=True)
        return _with_annplyr_accessor(
            relocate_adata(
                self._obj,
                **kwargs,
                before=relocate_anchor(self._obj, before, kwargs),
                after=relocate_anchor(self._obj, after, kwargs),
                inplace=inplace,
            )
        )

    def distinct(
        self,
        *names: Selector,
        keep_all: bool = False,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        """Keep the first row for each combination of the named columns.

        The names come from one container. ``var`` names make features the rows.
        """
        kwargs = container_kwargs(self._obj, names, verb="distinct", allowed=META_AND_X, negation=True)
        if len(kwargs) > 1:
            raise IncompatibleAxisError("distinct takes columns from one container at a time")
        return _with_annplyr_accessor(
            distinct_adata(
                self._obj,
                **kwargs,
                axis="var" if "var" in kwargs else "obs",
                keep_all=keep_all,
                copy=copy,
                max_matrix_values=max_matrix_values,
            )
        )

    def left_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        relationship: JoinRelationship = "many-to-one",
        multiple: JoinMultiple = "error",
        unmatched: JoinUnmatched = "drop",
        na_matches: NaMatches = "na",
        suffixes: tuple[str, str] = ("", "_right"),
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            left_join_adata(
                self._obj,
                other,
                by=by,
                axis=axis,
                relationship=relationship,
                multiple=multiple,
                unmatched=unmatched,
                na_matches=na_matches,
                suffixes=suffixes,
                copy=copy,
            )
        )

    def inner_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        relationship: JoinRelationship = "many-to-one",
        multiple: JoinMultiple = "error",
        unmatched: JoinUnmatched = "drop",
        na_matches: NaMatches = "na",
        suffixes: tuple[str, str] = ("", "_right"),
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            inner_join_adata(
                self._obj,
                other,
                by=by,
                axis=axis,
                relationship=relationship,
                multiple=multiple,
                unmatched=unmatched,
                na_matches=na_matches,
                suffixes=suffixes,
                copy=copy,
            )
        )

    def right_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        relationship: JoinRelationship = "many-to-one",
        multiple: JoinMultiple = "error",
        unmatched: JoinUnmatched = "error",
        na_matches: NaMatches = "na",
        suffixes: tuple[str, str] = ("", "_right"),
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            right_join_adata(
                self._obj,
                other,
                by=by,
                axis=axis,
                relationship=relationship,
                multiple=multiple,
                unmatched=unmatched,
                na_matches=na_matches,
                suffixes=suffixes,
                copy=copy,
            )
        )

    def full_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        relationship: JoinRelationship = "many-to-one",
        multiple: JoinMultiple = "error",
        unmatched: JoinUnmatched = "error",
        na_matches: NaMatches = "na",
        suffixes: tuple[str, str] = ("", "_right"),
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            full_join_adata(
                self._obj,
                other,
                by=by,
                axis=axis,
                relationship=relationship,
                multiple=multiple,
                unmatched=unmatched,
                na_matches=na_matches,
                suffixes=suffixes,
                copy=copy,
            )
        )

    def semi_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        na_matches: NaMatches = "na",
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            semi_join_adata(self._obj, other, by=by, axis=axis, na_matches=na_matches, copy=copy)
        )

    def anti_join(
        self,
        other: JoinInput,
        *,
        by: JoinBy = None,
        axis: Axis = "obs",
        na_matches: NaMatches = "na",
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            anti_join_adata(self._obj, other, by=by, axis=axis, na_matches=na_matches, copy=copy)
        )

    def arrange(
        self,
        *by: Expression | Sequence[Expression],
        source: str | None = None,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            arrange_adata(
                self._obj,
                *by,
                source=source,
                copy=copy,
                max_matrix_values=max_matrix_values,
            )
        )

    def slice(self, *indices: Any, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_adata(self._obj, *indices, axis=axis, copy=copy))

    def slice_head(self, n: int = 5, *, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_head_adata(self._obj, n=n, axis=axis, copy=copy))

    def slice_tail(self, n: int = 5, *, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_tail_adata(self._obj, n=n, axis=axis, copy=copy))

    def slice_min(self, by: Expression, n: int = 5, *, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_min_adata(self._obj, by=by, n=n, axis=axis, copy=copy))

    def slice_max(self, by: Expression, n: int = 5, *, axis: Axis = "obs", copy: bool = True) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(slice_max_adata(self._obj, by=by, n=n, axis=axis, copy=copy))

    def slice_sample(
        self,
        n: int | None = None,
        *,
        prop: float | None = None,
        replace: bool = False,
        random_state: int | None = None,
        axis: Axis = "obs",
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            slice_sample_adata(
                self._obj,
                n=n,
                prop=prop,
                replace=replace,
                random_state=random_state,
                axis=axis,
                copy=copy,
            )
        )

    def mutate(
        self,
        assignments: Mapping[str, Any] | None = None,
        /,
        *,
        into: str = "obs",
        inplace: bool = False,
        max_matrix_values: int | None = None,
        **named: Any,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            mutate_adata(
                self._obj,
                merge_named(assignments, named, verb="mutate"),
                into=into,
                inplace=inplace,
                max_matrix_values=max_matrix_values,
            )
        )

    def transmute(
        self,
        assignments: Mapping[str, Any] | None = None,
        /,
        *,
        into: str = "obs",
        max_matrix_values: int | None = None,
        **named: Any,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            transmute_adata(
                self._obj,
                merge_named(assignments, named, verb="transmute"),
                into=into,
                max_matrix_values=max_matrix_values,
            )
        )

    @overload
    def group_by(self) -> AnnDataWithAnnplyr: ...

    @overload
    def group_by(self, *by: Selector) -> GroupedAnnData: ...

    def group_by(self, *by: Selector) -> AnnDataWithAnnplyr | GroupedAnnData:
        """Group rows so later verbs run group-locally, one key per argument.

        The keys are positional, matching :meth:`count`, so grouping reads the
        same as counting::

            adata.ap.group_by(obs.sample, obs.cell_type).summarize({"n": ap.n()})

        A key is a string, a column such as ``obs.louvain`` or a selector.
        Computed values are not keys: ``mutate`` them first. ``var`` keys group
        features instead of observations, and mixing the two raises
        ``IncompatibleAxisError``. With no keys the object comes back unchanged,
        so a chain is never broken by a conditional ``group_by``.

        Parameters
        ----------
        *by
            See the annotated public signature.

        Returns
        -------
        AnnData or GroupedAnnData

        Ownership
        ---------
        The wrapper references the exact input object and never copies it.

        Raises
        ------
        SelectionError
            If the request violates the corresponding public contract.
        UnknownColumnError
            If the request violates the corresponding public contract.
        UnknownSourceError
            If the request violates the corresponding public contract.
        IncompatibleAxisError
            If the request violates the corresponding public contract.
        AnnplyrError
            If the request violates the corresponding public contract.
        """
        if not by:
            return _with_annplyr_accessor(self._obj)
        axis, keys = key_names(self._obj, by, verb="group_by")
        return GroupedAnnData(self._obj, obs=keys) if axis == "obs" else GroupedAnnData(self._obj, var=keys)

    def summarize(
        self,
        stats: Mapping[str, Any] | None = None,
        /,
        *,
        by: Selector | None = None,
        source: str | None = None,
        max_matrix_values: int | None = None,
        **named: Any,
    ) -> pd.DataFrame:
        return summarize_adata(
            self._obj,
            merge_named(stats, named, verb="summarize"),
            by=by,
            source=source,
            max_matrix_values=max_matrix_values,
        )

    summarise = summarize

    def count(
        self,
        *by: Selector,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis | None = None,
        name: str = "n",
    ) -> pd.DataFrame:
        resolved, keys = key_names(self._obj, by, verb="count", axis=axis, qualified=True)
        return count_adata(self._obj, by=keys or None, wt=wt, sort=sort, axis=resolved, name=name)

    def tally(
        self,
        *by: Selector,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis | None = None,
        name: str = "n",
    ) -> pd.DataFrame:
        resolved, keys = key_names(self._obj, by, verb="tally", axis=axis, qualified=True)
        return tally_adata(self._obj, by=keys or None, wt=wt, sort=sort, axis=resolved, name=name)

    def add_count(
        self,
        *by: Selector,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis | None = None,
        name: str = "n",
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        resolved, keys = key_names(self._obj, by, verb="add_count", axis=axis, qualified=True)
        return _with_annplyr_accessor(
            add_count_adata(self._obj, by=keys or None, wt=wt, sort=sort, axis=resolved, name=name, inplace=inplace)
        )

    def add_tally(
        self,
        *,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis = "obs",
        name: str = "n",
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            add_tally_adata(self._obj, wt=wt, sort=sort, axis=axis, name=name, inplace=inplace)
        )

    def pull(self, name: Selector, *, max_matrix_values: int | None = None) -> pd.Series:
        """Return one column as a Series, such as ``pull(obs.louvain)`` or ``pull(counts.MS4A1)``."""
        kwargs = container_kwargs(self._obj, [name], verb="pull", allowed=_PULLABLE, whole_matrix=False)
        return pull_adata(self._obj, **kwargs, max_matrix_values=max_matrix_values)

    def to_wide(
        self,
        *names: Selector,
        allow_all_features: bool = False,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        """Return the named columns from any container as one DataFrame, one row per cell.

        A whole matrix (``X``, ``raw`` or a layer) needs ``allow_all_features=True``.
        A whole ``obs`` or ``obsm`` key does not.
        """
        kwargs = container_kwargs(self._obj, names, verb="to_wide", allowed=_TABLE, whole_matrix=allow_all_features)
        return to_wide_adata(self._obj, **kwargs, max_matrix_values=max_matrix_values)

    def glimpse(self, n: int = 5) -> None:
        """Print the first rows of ``obs`` and the names of layers, obsm keys and raw."""
        adata = self._obj
        listed = lambda keys: ", ".join(str(k) for k in keys if k is not None) or "none"
        print(f"{adata.n_obs} cells x {adata.n_vars} features")
        print(cast(pd.DataFrame, adata.obs).head(n).to_string())
        print(f"layers: {listed(adata.layers.keys())}")
        print(f"obsm: {listed(adata.obsm.keys())}")
        print(f"raw: {'yes' if adata.raw is not None else 'no'}")

    def to_tidy(
        self,
        obs: Selector | None = None,
        x: Selector | None = None,
        raw: Selector | None = None,
        *,
        layer: str | None = None,
        obs_name: str = "obs_name",
        feature: str = "feature",
        value: str = "value",
        allow_all_features: bool = False,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return to_tidy_adata(
            self._obj,
            obs=obs,
            x=x,
            raw=raw,
            layer=layer,
            obs_name=obs_name,
            feature=feature,
            value=value,
            allow_all_features=allow_all_features,
            max_matrix_values=max_matrix_values,
        )

    def pivot_longer(
        self,
        *names: Selector,
        obs_name: str = "obs_name",
        names_to: str = "name",
        values_to: str = "value",
        allow_all_features: bool = False,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        """Return one row per cell and feature: ``obs`` names ride along, matrix names are pivoted.

        A whole matrix (``X``, ``raw`` or a layer) needs ``allow_all_features=True``.
        """
        kwargs = container_kwargs(
            self._obj,
            names,
            verb="pivot_longer",
            allowed=frozenset({"obs", "x", "raw", "layer"}),
            whole_matrix=allow_all_features,
        )
        return pivot_longer_adata(
            self._obj,
            **kwargs,
            obs_name=obs_name,
            names_to=names_to,
            values_to=values_to,
            allow_all_features=allow_all_features,
            max_matrix_values=max_matrix_values,
        )

    def as_frame(
        self,
        source: Source,
        *,
        key: str | None = None,
        select: Selector | None = None,
        layer: str | None = None,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return as_frame_adata(
            self._obj, source=source, key=key, select=select, layer=layer, max_matrix_values=max_matrix_values
        )

    def nest_by(
        self,
        *,
        by: Selector,
        obs: Selector | None = None,
        var: Selector | None = None,
        axis: Axis = "obs",
        name: str = "data",
    ) -> pd.DataFrame:
        return nest_by_adata(self._obj, by=by, obs=obs, var=var, axis=axis, name=name)

    def pipe(self, func: Callable[..., T] | tuple[Callable[..., T], str], *args: Any, **kwargs: Any) -> T:
        if isinstance(func, tuple):
            call, data_keyword = func
            if data_keyword in kwargs:
                msg = f"{data_keyword!r} is both the pipe target and a keyword argument"
                raise AnnplyrError(msg)
            kwargs[data_keyword] = self._obj
            return call(*args, **kwargs)
        return func(self._obj, *args, **kwargs)


install_verbs(AnnplyrAccessor)
