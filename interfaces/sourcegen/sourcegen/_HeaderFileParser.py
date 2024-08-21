"""Parser for existing CLib headers."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging
import re
from typing import List

from ._dataclasses import HeaderFile, Func


_logger = logging.getLogger()

_clib_path = Path(__file__).parent.joinpath("../../../include/cantera/clib").resolve()
_clib_ignore = ["clib_defs.h", "ctmatlab.h"]

class HeaderFileParser:

    def __init__(self, path: Path, ignore_funcs: List[str] = None):
        self._path = path
        self._ignore_funcs = ignore_funcs

    @classmethod
    def from_headers(cls, ignore_files, ignore_funcs) -> List[HeaderFile]:
        """Parse existing header file."""
        files = [_ for _ in _clib_path.glob("*.h")
                 if _.name not in ignore_files + _clib_ignore]
        return [cls(_, ignore_funcs.get(_.name, []))._parse_h() for _ in files]

    def _parse_h(self) -> HeaderFile:
        ct = self._path.read_text()

        matches = re.finditer(r"CANTERA_CAPI.*?;", ct, re.DOTALL)
        c_functions = [re.sub(r"\s+", " ", m.group()).replace("CANTERA_CAPI ", "")
                       for m in matches]

        if not c_functions:
            return

        parsed = map(Func.from_str, c_functions)

        _logger.info(f"  parsing {self._path.name!r}")
        if self._ignore_funcs:
            _logger.info(f"    ignoring {self._ignore_funcs!r}")

        parsed = [f for f in parsed if f.name not in self._ignore_funcs]

        if not parsed:
            return

        return HeaderFile(self._path, parsed)
