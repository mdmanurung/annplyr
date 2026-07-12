# Adapted from the scanpy documentation extension for compact typed returns.
from __future__ import annotations

import re
from collections.abc import Generator, Iterable

from sphinx.application import Sphinx
from sphinx.ext.napoleon import NumpyDocstring


def _process_return(lines: Iterable[str]) -> Generator[str]:
    for line in lines:
        if match := re.fullmatch(r"(?P<param>\w+)\s+:\s+(?P<type>[\w.]+)", line):
            yield f"-{match['param']} (:class:`~{match['type']}`)"
        else:
            yield line


def _parse_returns_section(self: NumpyDocstring, section: str) -> list[str]:
    lines_raw = self._dedent(self._consume_to_next_section())
    if lines_raw and lines_raw[0] == ":":
        del lines_raw[0]
    lines = self._format_block(":returns: ", list(_process_return(lines_raw)))
    if lines and lines[-1]:
        lines.append("")
    return lines


def setup(app: Sphinx) -> None:
    """Patch napoleon's NumPy returns parser."""
    NumpyDocstring._parse_returns_section = _parse_returns_section
