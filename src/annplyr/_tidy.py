"""Chainable AnnData result type shared by the accessor and function surfaces.

Every AnnData-returning verb hands its result through :func:`TidyAnnData.wrap`,
so the result keeps the verbs and a pipeline never has to re-enter ``.ap``.
``TidyAnnData`` subclasses :class:`anndata.AnnData`, so ``isinstance`` checks in
Scanpy, AnnData, and user code keep passing.
"""

from __future__ import annotations

import inspect
from collections.abc import Callable
from typing import TYPE_CHECKING, Any

from anndata import AnnData

#: Composed ``TidyAnnData + user subclass`` types, cached per user class so the
#: dynamic class is created once. Also registered in this module's namespace so
#: in-process pickling can resolve the name.
_COMPOSED: dict[type, type] = {}

#: Maps a composed class back to the user class it was built from, so
#: :meth:`TidyAnnData.to_anndata` can restore it.
_ORIGINAL_CLASS: dict[type, type] = {}

if TYPE_CHECKING:
    from annplyr._accessor import AnnplyrAccessor


def _restore(result: Any) -> Any:
    """Undo the verb upgrade, returning the original AnnData class.

    A plain ``TidyAnnData`` goes back to ``AnnData``; a composed
    ``TidyAnnData + user subclass`` goes back to the user's class.
    """
    if not isinstance(result, AnnData):
        return result
    if type(result) is TidyAnnData:
        original: type = AnnData
    else:
        composed_original = _ORIGINAL_CLASS.get(type(result))
        if composed_original is None:
            return result
        original = composed_original
    try:
        result.__class__ = original
    except TypeError:
        pass
    return result


def _reclassify(source: Any, result: Any) -> Any:
    """Carry ``source``'s user subclass onto ``result`` and re-add the verbs.

    ``AnnData`` builds fresh objects as plain ``AnnData`` inside ``copy``,
    ``__getitem__``, and friends, so a user subclass has to be reapplied
    explicitly before the verbs are added back.
    """
    if not isinstance(result, AnnData):
        return result
    original = _ORIGINAL_CLASS.get(type(source))
    if original is not None:
        try:
            result.__class__ = original
        except TypeError:
            return result
    return TidyAnnData.wrap(result)


def _upgrade_target(current: type) -> type | None:
    """Return the class ``current`` should be upgraded to, or ``None``.

    ``None`` means the object cannot be upgraded without losing information, in
    which case callers must leave it alone.
    """
    if current is AnnData:
        return TidyAnnData
    cached = _COMPOSED.get(current)
    if cached is not None:
        return cached
    try:
        composed = type(f"Tidy{current.__name__}", (TidyAnnData, current), {})
    except TypeError:
        return None
    _COMPOSED[current] = composed
    _ORIGINAL_CLASS[composed] = current
    globals()[composed.__name__] = composed
    return composed


class TidyAnnData(AnnData):
    """An ``AnnData`` that also carries the annplyr verbs.

    Use :meth:`wrap` rather than calling the constructor. It upgrades an
    existing object in place, so the object keeps its identity.

    A user-defined ``AnnData`` subclass keeps its own methods, slots, and
    ``isinstance`` relationships: the upgrade composes a class with it rather
    than replacing it. Verbs that derive a *new* object (``filter``, ``select``, and friends) return a plain ``TidyAnnData``, which
    matches upstream ``AnnData``, whose own ``copy`` and ``__getitem__`` also
    return a plain ``AnnData`` for subclasses.

    The static declaration of every verb lives in the generated
    ``annplyr/_tidy.pyi`` stub, not here, so there is one source of truth.
    """

    @classmethod
    def wrap(cls, source: Any) -> Any:
        """Return ``source`` carrying the annplyr verbs, preserving identity.

        Non-AnnData values (the ``pandas`` results of terminal verbs) and
        already-wrapped objects are returned untouched. When ``source`` is an
        instance of a user-defined ``AnnData`` subclass, the upgrade composes a
        class with that subclass so the user's methods, slots, and
        ``isinstance`` relationships survive. If composition or the in-place
        class change is not possible, ``source`` is returned unchanged rather
        than downgraded or raising.
        """
        if isinstance(source, cls) or not isinstance(source, AnnData):
            return source
        target = _upgrade_target(type(source))
        if target is None:
            return source
        try:
            source.__class__ = target
        except TypeError:
            return source
        return source

    def to_anndata(self) -> AnnData:
        """Return an independent ``AnnData`` without the annplyr verbs.

        A user-defined ``AnnData`` subclass this object still carries is
        restored rather than erased, so removing the verbs never silently drops
        a subclass the object already had.
        """
        return _restore(self.copy())

    def write_zarr(self, *args: Any, **kwargs: Any) -> None:
        # anndata's IO registry keys on the exact class and has no zarr writer for a
        # subclass, so write as the class the object had before the verbs were added.
        current = type(self)
        plain = _ORIGINAL_CLASS.get(current, AnnData)
        self.__class__ = plain
        try:
            plain.write_zarr(self, *args, **kwargs)
        finally:
            self.__class__ = current

    def __getitem__(self, index: Any) -> Any:
        return _reclassify(self, super().__getitem__(index))

    def copy(self, *args: Any, **kwargs: Any) -> Any:
        return _reclassify(self, super().copy(*args, **kwargs))

    def to_memory(self, *args: Any, **kwargs: Any) -> Any:
        return _reclassify(self, super().to_memory(*args, **kwargs))

    def rename_categories(self, *args: Any, **kwargs: Any) -> Any:
        return _reclassify(self, super().rename_categories(*args, **kwargs))

    if hasattr(AnnData, "concatenate"):  # removed in anndata 0.13

        def concatenate(self, *args: Any, **kwargs: Any) -> Any:
            return _reclassify(self, super().concatenate(*args, **kwargs))


def _without_self(signature: inspect.Signature) -> inspect.Signature:
    return signature.replace(parameters=[p for p in signature.parameters.values() if p.name != "self"])


def _delegate_to_accessor(name: str, accessor_method: Callable[..., Any]) -> Callable[..., Any]:
    def verb(self: TidyAnnData, *args: Any, **kwargs: Any) -> Any:
        return getattr(self.ap, name)(*args, **kwargs)

    verb.__name__ = name
    verb.__qualname__ = f"TidyAnnData.{name}"
    verb.__module__ = __name__
    verb.__doc__ = accessor_method.__doc__
    verb.__signature__ = _without_self(inspect.signature(accessor_method))  # type: ignore[attr-defined]
    return verb


#: Accessor verbs that could not be attached because AnnData already owns the
#: name. Kept so tests fail loudly if a future anndata release adds a collision
#: rather than silently dropping a verb.
SHADOWED_VERBS: dict[str, str] = {}

#: Names already attached, so a repeated :func:`install_verbs` call is a no-op
#: instead of reporting its own verbs as shadowed.
_INSTALLED: set[str] = set()


def install_verbs(accessor_cls: type[AnnplyrAccessor]) -> None:
    """Attach every public accessor verb to :class:`TidyAnnData` as a method.

    Verbs are delegated rather than redeclared so the accessor stays the single
    source of truth for signatures, docstrings, and behavior. Calling this more
    than once is safe and does not re-report already-attached verbs.
    """
    for name, member in inspect.getmembers(accessor_cls, predicate=inspect.isfunction):
        if name.startswith("_"):
            continue
        if name in _INSTALLED:
            continue
        if hasattr(TidyAnnData, name):
            SHADOWED_VERBS[name] = "already defined on AnnData"
            continue
        setattr(TidyAnnData, name, _delegate_to_accessor(name, member))
        _INSTALLED.add(name)
