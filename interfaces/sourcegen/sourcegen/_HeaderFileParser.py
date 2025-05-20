"""Parser for YAML header configurations or existing CLib headers."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging
import re
from sys import version_info

if version_info.minor < 11:
    from typing import Iterable
    from typing_extensions import Self
else:
    from typing import Iterable, Self

from ._dataclasses import HeaderFile, Func, Recipe
from ._helpers import read_config
from .headers import RecipeParser


_LOGGER = logging.getLogger()

_CLIB_PATH = Path(__file__).parents[3] / "include" / "cantera" / "clib"
_CLIB_IGNORE = ["clib_defs.h", "ctmatlab.h"]

_HERE = Path(__file__).parent

class HeaderFileParser:
    """
    Parser for header files or corresponding YAML specifications.

    Provides for convenience methods to generate lists of `HeaderFile` objects, which
    themselves are used for subsequent code scaffolding.
    """

    def __init__(self, path: Path, ignore_funcs: Iterable[str] = None) -> None:
        self._path = path
        self._ignore_funcs = ignore_funcs

    @classmethod
    def headers_from_yaml(
            cls: Self, ignore_files: Iterable[str], ignore_funcs: Iterable[str]
        ) -> list[HeaderFile]:
        """Parse header file YAML configuration."""
        files = sorted(
            ff for ff in (_HERE / "headers").glob("ct*.yaml")
            if ff.name not in ignore_files)
        files = [cls(ff, ignore_funcs.get(ff.name, []))._parse_yaml() for ff in files]
        bases = [f.base for f in files]

        # preprocess header information (uses CLibSourceGenerator)
        config = read_config(_HERE / "headers" / "config.yaml")
        templates = read_config(_HERE / "headers" / "templates.yaml")
        for key in ["ignore_files", "ignore_funcs"]:
            config.pop(key)
        RecipeParser(config, templates, bases).resolve_tags(files)
        # CLibSourceGenerator(None, config, templates).resolve_tags(files)
        return files

    def _parse_yaml(self) -> HeaderFile:
        msg = f"  parsing {self._path.name!r}"
        _LOGGER.info(msg)
        config = read_config(self._path)
        if self._ignore_funcs:
            msg = f"    ignoring {self._ignore_funcs!r}"
            _LOGGER.info(msg)

        recipes = []
        docstring = config["docstring"].split("\n")
        prefix = config["prefix"]
        base = config["base"]
        parents = config.get("parents", [])
        derived = config.get("derived", {})
        for recipe in config["recipes"]:
            if recipe["name"] in self._ignore_funcs:
                continue
            uses = recipe.get("uses", [])
            if not isinstance(uses, list):
                uses = [uses]
            recipes.append(
                Recipe(recipe["name"],
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
                          docstring)

    @classmethod
    def headers_from_h(
            cls: Self, ignore_files: Iterable[str], ignore_funcs: Iterable[str]
        ) -> list[HeaderFile]:
        """Parse existing header file."""
        files = [ff for ff in _CLIB_PATH.glob("*.h")
                 if ff.name not in ignore_files + _CLIB_IGNORE]
        files.sort()
        return [cls(ff, ignore_funcs.get(ff.name, []))._parse_h() for ff in files]

    def _parse_h(self) -> HeaderFile:
        ct = self._path.read_text(encoding="utf-8")

        matches = re.finditer(r"CANTERA_CAPI.*?;", ct, re.DOTALL)
        c_functions = [re.sub(r"\s+", " ", m.group()).replace("CANTERA_CAPI ", "")
                       for m in matches]

        if not c_functions:
            return

        parsed = map(Func.from_str, c_functions)

        msg = f"  parsing {self._path.name!r}"
        _LOGGER.info(msg)
        if self._ignore_funcs:
            msg = f"    ignoring {self._ignore_funcs!r}"
            _LOGGER.info(msg)

        parsed = [f for f in parsed if f.name not in self._ignore_funcs]

        if not parsed:
            return

        return HeaderFile(self._path, parsed)
