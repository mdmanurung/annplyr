"""Free-function forms of the annplyr verbs.

Every verb is available as a plain function taking the AnnData object first, so
a pipeline can start without the accessor and then chain with methods:

    import annplyr as ap

    summary = (
        ap.filter(adata, obs=ap.col("n_genes") >= 500)
        .select(obs=ap.everything(), x=panel)
        .summarize(obs={"cells": ap.n()}, by="louvain")
    )

The functions are generated from :class:`~annplyr.AnnplyrAccessor`, so there is
one implementation of each verb and the accessor, the class methods, and these
functions can never drift apart.
"""

from __future__ import annotations

import inspect
from collections.abc import Callable
from typing import Any, cast

from anndata import AnnData

from annplyr._accessor import AnnplyrAccessor


def _make_function(name: str, accessor_method: Callable[..., Any]) -> Callable[..., Any]:
    def verb(adata: AnnData, *args: Any, **kwargs: Any) -> Any:
        return getattr(cast(Any, adata).ap, name)(*args, **kwargs)

    verb.__name__ = name
    verb.__qualname__ = name
    verb.__module__ = __name__
    verb.__doc__ = (
        f"Apply ``{name}`` to an AnnData object.\n\n"
        f"The first argument is the AnnData object, matching the accessor call\n"
        f"with the object moved in front. Every other argument is identical to\n"
        f"``AnnData.ap.{name}``, and the result carries the verbs, so this call\n"
        f"chains with the following step.\n\n"
        f"{accessor_method.__doc__ or ''}"
    ).rstrip() + "\n"
    verb.__annotations__ = {"adata": AnnData, "return": Any}
    parameters = [inspect.Parameter("adata", inspect.Parameter.POSITIONAL_OR_KEYWORD, annotation=AnnData)]
    parameters += [p for p in inspect.signature(accessor_method).parameters.values() if p.name != "self"]
    verb.__signature__ = inspect.Signature(parameters)  # type: ignore[attr-defined]
    return verb


#: Free-function verb surface, keyed by verb name.
FUNCTIONS: dict[str, Callable[..., Any]] = {
    name: _make_function(name, member)
    for name, member in inspect.getmembers(AnnplyrAccessor, predicate=inspect.isfunction)
    if not name.startswith("_")
}

globals().update(FUNCTIONS)

__all__ = sorted(FUNCTIONS)
