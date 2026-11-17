"""Plot-ready verbs must hand plotnine a frame it can draw.

annplyr has no runtime plotnine dependency, but its docs and ggann draw its
tables with plotnine. plotnine cannot take pandas sparse columns, so every verb
that the docs call "plot-ready" has to return dense columns, even when ``X`` is
sparse. This file is the data contract between annplyr and the plot packages.

It skips when plotnine is absent. annplyr's own test environment does not
install plotnine, so it runs only in the downstream job of plotnine-extra's CI,
or after installing the ``doc`` group.
"""

from __future__ import annotations

import pandas as pd
import pytest

import annplyr as ap
from annplyr.data import X, obs

p9 = pytest.importorskip("plotnine")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")


def _to_wide(adata):
    return adata.ap.to_wide(obs.batch, obs.score, X.g0, X.g2)


def _join_features(adata):
    return adata.ap.join_features(X.g0)


def _summarize(adata):
    return adata.ap.summarize(cells=ap.n(), mean_score=ap.mean(obs.score), by=obs.batch)


def _count(adata):
    return adata.ap.count("batch")


CASES = {
    "to_wide": _to_wide,
    "join_features": _join_features,
    "summarize": _summarize,
    "count": _count,
}

VERBS = list(CASES)


def _xy(frame: pd.DataFrame) -> tuple[str, str]:
    numeric = [c for c in frame.columns if pd.api.types.is_numeric_dtype(frame[c])]
    if len(numeric) >= 2:
        return numeric[0], numeric[1]
    return str(frame.columns[0]), numeric[0]


@pytest.mark.parametrize("verb", VERBS)
def test_output_has_no_sparse_columns(verb, sparse_adata):
    frame = CASES[verb](sparse_adata)
    sparse = [c for c in frame.columns if isinstance(frame[c].dtype, pd.SparseDtype)]
    assert sparse == [], f"{verb} returned sparse columns plotnine cannot draw: {sparse}"


@pytest.mark.parametrize("verb", VERBS)
def test_output_draws_with_plotnine(verb, sparse_adata):
    frame = CASES[verb](sparse_adata)
    x, y = _xy(frame)
    plot = p9.ggplot(frame, p9.aes(x, y)) + p9.geom_point()
    p9.ggplot.draw(plot, show=False)
