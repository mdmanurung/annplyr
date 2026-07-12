from __future__ import annotations

import warnings

from anndata import register_anndata_namespace


class AccessorRegistrationWarning(UserWarning):
    """Warning issued when a conflicting accessor is registered."""


def register_anndata_accessor(name: str):
    """Register an AnnData namespace using AnnData's public extension API."""

    def decorator(accessor: type):
        try:
            return register_anndata_namespace(name)(accessor)
        except AttributeError:
            warnings.warn(
                f"registration of accessor {accessor!r} under name {name!r} is overriding a preexisting attribute",
                AccessorRegistrationWarning,
                stacklevel=2,
            )
            return accessor

    return decorator
