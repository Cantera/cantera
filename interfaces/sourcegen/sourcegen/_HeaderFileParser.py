"""Parser for YAML header configurations or existing CLib headers."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging
import re

from ._dataclasses import HeaderFile, Func, Recipe
from ._helpers import read_config


_LOGGER = logging.getLogger()

_CLIB_PATH = Path(__file__).parents[3] / "include" / "cantera" / "clib"
_CLIB_IGNORE = ["clib_defs.h", "ctmatlab.h"]

_DATA_PATH = Path(__file__).parent / "_data"

class HeaderFileParser:
    """
    Parser for header files or corresponding YAML specifications.

    Provides for convenience methods to generate lists of `HeaderFile` objects, which
    themselves are used for subsequent code scaffolding.
    """

    def __init__(self, path: Path, ignore_funcs: list[str] = None):
        self._path = path
        self._ignore_funcs = ignore_funcs

    @classmethod
    def headers_from_yaml(cls, ignore_files, ignore_funcs) -> list[HeaderFile]:
        """Parse header file YAML configuration."""
        files = [ff for ff in _DATA_PATH.glob("*.yaml") if ff.name not in ignore_files]
        files.sort()
        return [cls(ff, ignore_funcs.get(ff.name, []))._parse_yaml() for ff in files]

    def _parse_yaml(self) -> HeaderFile:
        def read_docstring():
            doc = []
            with self._path.open("r", encoding="utf-8") as fid:
                while True:
                    line = fid.readline()
                    if line.startswith("#"):
                        doc.append(line.lstrip("#").strip())
                    else:
                        break
            if doc and doc[0].startswith("This file is part of "):
                return []
            return doc

        config = read_config(self._path)
        recipes = []
        prefix = config["prefix"]
        base = config["base"]
        parents = config.get("parents", [])
        derived = config.get("derived", [])
        for recipe in config["recipes"]:
            if recipe['name'] in self._ignore_funcs:
                continue
            uses = recipe.get("uses", [])
            if not isinstance(uses, list):
                uses = [uses]
            recipes.append(
                Recipe(recipe['name'],
                       recipe.get("implements", ""),
                       uses,
                       recipe.get("what", ""),
                       recipe.get("brief", ""),
                       recipe.get("code", ""),
                       prefix,
                       base,
                       parents,
                       derived))

        return HeaderFile(self._path, [], prefix, base, parents, derived, recipes,
                          read_docstring())

    @classmethod
    def headers_from_h(cls, ignore_files, ignore_funcs) -> list[HeaderFile]:
        """Parse existing header file."""
        files = [ff for ff in _CLIB_PATH.glob("*.h")
                 if ff.name not in ignore_files + _CLIB_IGNORE]
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

        _LOGGER.info(f"  parsing {self._path.name!r}")
        if self._ignore_funcs:
            _LOGGER.info(f"    ignoring {self._ignore_funcs!r}")

        parsed = [f for f in parsed if f.name not in self._ignore_funcs]

        if not parsed:
            return

        return HeaderFile(self._path, parsed)
