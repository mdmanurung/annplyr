from __future__ import annotations


class AnnplyrError(Exception):
    """Base class for annplyr contract errors."""


class SelectionError(AnnplyrError, ValueError):
    """Raised when a selector is invalid or resolves to an invalid shape."""


class UnknownColumnError(SelectionError, KeyError):
    """Raised when a selector references a missing column."""


class UnknownSourceError(AnnplyrError, KeyError):
    """Raised when an AnnData source, layer, or aligned mapping key is missing."""


class DuplicateNameError(SelectionError):
    """Raised when duplicate names make an operation ambiguous."""


class NameRepairError(SelectionError):
    """Raised when annplyr cannot repair generated names safely."""


class IncompatibleAxisError(AnnplyrError, ValueError):
    """Raised when a verb receives incompatible obs/var axis inputs."""


class SizeMismatchError(AnnplyrError, ValueError):
    """Raised when computed values cannot be aligned to an AnnData axis."""


class JoinRelationshipError(AnnplyrError, ValueError):
    """Raised when a join violates the requested relationship contract."""
