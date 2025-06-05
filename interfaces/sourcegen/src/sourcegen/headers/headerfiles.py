"""Parser for YAML header configurations."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging
from sys import version_info

if version_info.minor < 11:
    from typing import Iterable
    from typing_extensions import Self
else:
    from typing import Iterable, Self

from ..dataclasses import HeaderFile, Recipe
from .._helpers import read_config
from .generator import HeaderGenerator


_LOGGER = logging.getLogger()

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
            cls: Self, ignore_files: Iterable[str], ignore_funcs: Iterable[str],
            root: str = ".",
        ) -> list[HeaderFile]:
        """Parse header file YAML configuration."""
        files = sorted(
            ff for ff in _HERE.glob("ct*.yaml")
            if ff.name not in ignore_files)
        files = [cls(ff, ignore_funcs.get(ff.name, []))._parse_yaml() for ff in files]
        bases = [f.base for f in files]

        # preprocess header information (uses CLibSourceGenerator)
        config = read_config(_HERE / "config.yaml")
        templates = read_config(_HERE / "templates.yaml")
        for key in ["ignore_files", "ignore_funcs"]:
            config.pop(key)
        HeaderGenerator(config, templates, bases).resolve_tags(files, root)
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
                       recipe.get("brief", ""),
                       recipe.get("what", ""),
                       recipe.get("declaration", ""),
                       recipe.get("parameters", ""),
                       recipe.get("returns", ""),
                       uses,
                       recipe.get("wraps", ""),
                       recipe.get("code", ""),
                       prefix,
                       base,
                       parents,
                       derived))

        return HeaderFile(self._path, [], prefix, base, parents, derived, recipes,
                          docstring)
