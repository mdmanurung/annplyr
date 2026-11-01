from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any, TypeVar, cast, overload

import pandas as pd
from anndata import AnnData

from annplyr._errors import AnnplyrError
from annplyr._extensions import register_anndata_accessor
from annplyr._grouped import GroupedAnnData
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
    SourceSelectors,
)
from annplyr._verbs import (
    _axis,
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

    def select(
        self,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        copy: bool = True,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(select_adata(self._obj, obs=obs, var=var, x=x, copy=copy))

    def rename(
        self,
        obs: Mapping[str, str] | None = None,
        var: Mapping[str, str] | None = None,
        x: Mapping[str, str] | None = None,
        *,
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(rename_adata(self._obj, obs=obs, var=var, x=x, inplace=inplace))

    def rename_with(
        self,
        func: Callable[[str], str],
        *,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(rename_with_adata(self._obj, func, obs=obs, var=var, x=x, inplace=inplace))

    def relocate(
        self,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        *,
        before: str | None = None,
        after: str | None = None,
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            relocate_adata(self._obj, obs=obs, var=var, x=x, before=before, after=after, inplace=inplace)
        )

    def distinct(
        self,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        *,
        axis: Axis = "obs",
        keep_all: bool = False,
        copy: bool = True,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            distinct_adata(
                self._obj,
                obs=obs,
                var=var,
                x=x,
                axis=axis,
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
        *,
        into: str = "obs",
        inplace: bool = False,
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            mutate_adata(
                self._obj,
                assignments,
                into=into,
                inplace=inplace,
                max_matrix_values=max_matrix_values,
            )
        )

    def transmute(
        self,
        assignments: Mapping[str, Any] | None = None,
        *,
        into: str = "obs",
        max_matrix_values: int | None = None,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            transmute_adata(
                self._obj,
                assignments,
                into=into,
                max_matrix_values=max_matrix_values,
            )
        )

    @overload
    def group_by(self) -> AnnDataWithAnnplyr: ...

    @overload
    def group_by(self, *by: Selector, axis: Axis = "obs") -> GroupedAnnData: ...

    def group_by(self, *by: Selector, axis: Axis = "obs") -> AnnDataWithAnnplyr | GroupedAnnData:
        """Group rows so later verbs run group-locally, one key per argument.

        The keys are positional, matching :meth:`count` and the ``by=`` of
        :meth:`summarize`, so grouping reads the same as summarizing::

            adata.ap.group_by("sample", "cell_type").summarize({"n": ap.n()})

        A bare name is ``obs`` metadata; qualify a feature name with ``var:``.
        With no keys the object comes back unchanged, so a chain is never broken
        by a conditional ``group_by``. ``axis="var"`` groups features instead of
        observations.

        Parameters
        ----------
        *by
            See the annotated public signature.
        axis
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
        axis = _axis(axis)
        if not by:
            return _with_annplyr_accessor(self._obj)
        if axis == "var":
            return GroupedAnnData(self._obj, var=cast(Selector, by))
        return GroupedAnnData(self._obj, obs=cast(Selector, by))

    def summarize(
        self,
        stats: Mapping[str, Any] | None = None,
        *,
        by: Selector | None = None,
        source: str | None = None,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return summarize_adata(
            self._obj,
            stats,
            by=by,
            source=source,
            max_matrix_values=max_matrix_values,
        )

    summarise = summarize

    def count(
        self,
        by: Selector | None = None,
        *,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis = "obs",
        name: str = "n",
    ) -> pd.DataFrame:
        return count_adata(self._obj, by=by, wt=wt, sort=sort, axis=axis, name=name)

    def tally(
        self,
        by: Selector | None = None,
        *,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis = "obs",
        name: str = "n",
    ) -> pd.DataFrame:
        return tally_adata(self._obj, by=by, wt=wt, sort=sort, axis=axis, name=name)

    def add_count(
        self,
        by: Selector | None = None,
        *,
        wt: Expression | None = None,
        sort: bool = False,
        axis: Axis = "obs",
        name: str = "n",
        inplace: bool = False,
    ) -> AnnDataWithAnnplyr:
        return _with_annplyr_accessor(
            add_count_adata(self._obj, by=by, wt=wt, sort=sort, axis=axis, name=name, inplace=inplace)
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

    def pull(
        self,
        obs: Selector | None = None,
        var: Selector | None = None,
        x: Selector | None = None,
        raw: Selector | None = None,
        obsm: SourceSelectors | None = None,
        varm: SourceSelectors | None = None,
        obsp: SourceSelectors | None = None,
        varp: SourceSelectors | None = None,
        uns: SourceSelectors | None = None,
        *,
        layer: str | None = None,
        max_matrix_values: int | None = None,
    ) -> pd.Series:
        return pull_adata(
            self._obj,
            obs=obs,
            var=var,
            x=x,
            raw=raw,
            obsm=obsm,
            varm=varm,
            obsp=obsp,
            varp=varp,
            uns=uns,
            layer=layer,
            max_matrix_values=max_matrix_values,
        )

    def to_wide(
        self,
        obs: Selector | None = None,
        x: Selector | None = None,
        raw: Selector | None = None,
        obsm: SourceSelectors | None = None,
        obsp: SourceSelectors | None = None,
        *,
        layer: str | None = None,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return to_wide_adata(
            self._obj,
            obs=obs,
            x=x,
            raw=raw,
            obsm=obsm,
            obsp=obsp,
            layer=layer,
            max_matrix_values=max_matrix_values,
        )

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
        obs: Selector | None = None,
        x: Selector | None = None,
        raw: Selector | None = None,
        *,
        layer: str | None = None,
        obs_name: str = "obs_name",
        names_to: str = "name",
        values_to: str = "value",
        allow_all_features: bool = False,
        max_matrix_values: int | None = None,
    ) -> pd.DataFrame:
        return pivot_longer_adata(
            self._obj,
            obs=obs,
            x=x,
            raw=raw,
            layer=layer,
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
