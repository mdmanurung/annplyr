"""Containers of an AnnData object as importable names.

``from annplyr.data import obs, var, X, raw`` gives the fixed containers. Any
other name is a layer or an ``obsm`` key, so ``from annplyr.data import counts``
works for ``adata.layers["counts"]``. A mistyped name fails at the verb with
``UnknownSourceError``, which lists the layers and keys that exist. A layer
named like an annplyr function, such as ``mean``, needs ``ap.layer("mean")``.
"""

from __future__ import annotations

from annplyr._expr import Source

obs = Source("obs")
var = Source("var")
X = Source("x")
raw = Source("raw")


def __getattr__(name: str) -> Source:
    import annplyr

    # Tools such as Sphinx probe modules for arbitrary names; refusing the
    # package name and the root exports keeps them from minting fake layers.
    if name.startswith("_") or name == "annplyr" or name in annplyr.__all__:
        msg = f"annplyr.data has no {name!r}; for a layer of that name use ap.layer({name!r})"
        raise AttributeError(msg)
    return Source(name)
