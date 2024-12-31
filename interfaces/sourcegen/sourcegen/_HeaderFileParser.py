"""Parser for YAML header configurations or existing CLib headers."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging
import re

from ._dataclasses import HeaderFile, Func, Recipe
from ._helpers import read_config

_logger = logging.getLogger()

_clib_path = Path(__file__).parents[3] / "include" / "cantera" / "clib"
_clib_ignore = ["clib_defs.h", "ctmatlab.h"]

_data_path = Path(__file__).parent / "_data"

class HeaderFileParser:

    def __init__(self, path: Path, ignore_funcs: list[str] = None):
        self._path = path
        self._ignore_funcs = ignore_funcs

    @classmethod
    def headers_from_yaml(cls, ignore_files, ignore_funcs) -> list[HeaderFile]:
        """Parse header file YAML configuration."""
        files = [ff for ff in _data_path.glob("*.yaml") if ff.name not in ignore_files]
        files.sort()
        return [cls(ff, ignore_funcs.get(ff.name, []))._parse_yaml() for ff in files]

    def _parse_yaml(self) -> HeaderFile:
        config = read_config(self._path)
        recipes = []
        cabinet = config.get("cabinet", [])
        prefix = cabinet["prefix"]
        base = cabinet["base"]
        parents = cabinet.get("parents", [])
        derived = cabinet.get("derived", [])
        uses = cabinet.get("uses", [])
        for func in cabinet["functions"]:
            if func['name'] in self._ignore_funcs:
                continue
            recipes.append(
                Recipe(prefix,
                       func['name'],
                       base,
                       parents,
                       derived,
                       uses,
                       func.get("implements", ""),
                       func.get("relates", []),
                       func.get("size-fcn", ""),
                       func.get("what", "")))
        return HeaderFile(self._path, [], recipes)

    @classmethod
    def headers_from_h(cls, ignore_files, ignore_funcs) -> list[HeaderFile]:
        """Parse existing header file."""
        files = [ff for ff in _clib_path.glob("*.h")
                 if ff.name not in ignore_files + _clib_ignore]
        files.sort()
        return [cls(ff, ignore_funcs.get(ff.name, []))._parse_h() for ff in files]

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
