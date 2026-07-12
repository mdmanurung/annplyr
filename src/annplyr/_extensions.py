from __future__ import annotations

import warnings

from anndata import AnnData


class AccessorRegistrationWarning(UserWarning):
    """Warning issued when a conflicting accessor is registered."""


class _CachedAccessor:
    def __init__(self, name: str, accessor: type):
        self._name = name
        self._accessor = accessor

    def __get__(self, obj, cls):
        if obj is None:
            return self._accessor

        try:
            cache = obj._cache
        except AttributeError:
            cache = obj._cache = {}

        try:
            return cache[self._name]
        except KeyError:
            pass

        accessor_obj = self._accessor(obj)
        cache[self._name] = accessor_obj
        return accessor_obj


def register_anndata_accessor(name: str):
    def decorator(accessor: type):
        if hasattr(AnnData, name):
            warnings.warn(
                f"registration of accessor {accessor!r} under name {name!r} is overriding a preexisting attribute",
                AccessorRegistrationWarning,
                stacklevel=2,
            )
        setattr(AnnData, name, _CachedAccessor(name, accessor))
        return accessor

    return decorator
