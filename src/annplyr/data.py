"""Containers of an AnnData object as importable names.

``from annplyr.data import obs, var, X, raw`` gives the fixed containers. Any
other name is a layer or an ``obsm`` key, so ``from annplyr.data import counts``
works for ``adata.layers["counts"]``. A mistyped name fails at the verb with
``UnknownSourceError``, which lists the layers and keys that exist.
"""

from __future__ import annotations

from annplyr._expr import Source

obs = Source("obs")
var = Source("var")
X = Source("x")
raw = Source("raw")


def __getattr__(name: str) -> Source:
    if name.startswith("_"):
        raise AttributeError(name)
    return Source(name)
