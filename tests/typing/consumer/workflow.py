from __future__ import annotations

from typing import assert_type, cast

import pandas as pd
from anndata import AnnData

import annplyr as ap
from annplyr.typing import (
    AnnDataWithAnnplyr,
    Expression,
    GroupedReturn,
    JoinInput,
    Selector,
    Source,
)


def extract_scores(adata: AnnDataWithAnnplyr) -> pd.DataFrame:
    return adata.ap.to_df(obs=["batch", "score"], x=["GATA3"])


def realistic_pipeline(adata: AnnData, metadata: pd.DataFrame) -> pd.DataFrame:
    typed = cast(AnnDataWithAnnplyr, adata)
    assert_type(typed.ap, ap.AnnplyrAccessor)

    predicate: Expression = ap.col("score") > 0
    columns: Selector = ap.all_of(["batch", "score"])
    source: Source = "x"
    join_input: JoinInput = metadata

    filtered = typed.ap.filter(obs=predicate)
    assert_type(filtered, AnnDataWithAnnplyr)
    selected = filtered.ap.select(obs=columns, x=["GATA3", "MS4A1"])
    assert_type(selected, AnnDataWithAnnplyr)

    expression = ap.col("score") + 1
    assert_type(expression, ap.AnnplyrExpr)
    mutated = selected.ap.mutate(obs={"score_plus_one": expression})
    assert_type(mutated, AnnDataWithAnnplyr)

    grouped_result: GroupedReturn = mutated.ap.group_by(obs="batch")
    assert_type(grouped_result, ap.GroupedAnnData | AnnDataWithAnnplyr)
    grouped = mutated.ap.group_by(obs="batch")
    assert_type(grouped, ap.GroupedAnnData)
    joined = grouped.left_join(join_input, by="batch")
    assert_type(joined, ap.GroupedAnnData)
    summary = joined.summarize(obs={"mean_score": ap.mean("score")})
    assert_type(summary, pd.DataFrame)

    ungrouped = joined.ungroup()
    assert_type(ungrouped, AnnDataWithAnnplyr)
    wide = ungrouped.ap.as_frame(source, select=["GATA3", "MS4A1"])
    assert_type(wide, pd.DataFrame)
    pulled = ungrouped.ap.pull(obs="score")
    assert_type(pulled, pd.Series)
    piped = ungrouped.ap.pipe(extract_scores)
    assert_type(piped, pd.DataFrame)
    return summary
