"""Parser for header YAML configurations or existing CLib headers."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging
import re
from typing import List

from ._dataclasses import HeaderFile, Func, Recipe, Param
from ._helpers import read_config


logger = logging.getLogger()

_clib_path = Path(__file__).parent.joinpath("../../../include/cantera/clib").resolve()
_clib_ignore = ["clib_defs.h", "ctmatlab.h"]
_data_path = Path(__file__).parent.joinpath("_data").resolve()

class HeaderFileParser:

    @staticmethod
    def _parse_param(c_param: str) -> Param:
        parts = c_param.strip().rsplit(" ", 1)
        if len(parts) == 2:
            return Param(*parts)

    @classmethod
    def _parse_func(cls, func_comment: tuple[str, str]) -> Func:
        c_func, annotations = func_comment
        lparen = c_func.index("(")
        rparen = c_func.index(")")
        front = c_func[0:lparen].split()

        params = list(map(cls._parse_param,
            filter(None, c_func[(lparen + 1):rparen].split(","))))

        ret_type = front[-2]
        name = front[-1]
        return Func(annotations, ret_type, name, params)

    def __init__(self, path: Path, ignore_funcs: List[str] = None):
        self._path = path
        self._ignore_funcs = ignore_funcs

    @classmethod
    def parse_yaml(cls, ignore_files, ignore_funcs) -> List[HeaderFile]:
        """Parse header file YAML configuration."""
        files = [_ for _ in _data_path.glob("*.yaml") if _.name not in ignore_files]
        return [cls(_, ignore_funcs.get(_.name, []))._parse_yaml() for _ in files]

    def _parse_yaml(self) -> HeaderFile:
        config = read_config(self._path)
        recipes = []
        for section in config.get("sections", []):
            prefix = section["prefix"]
            bases = section["bases"]
            for func in section["functions"]:
                func_name = f"{prefix}_{func['name']}"
                if func_name in self._ignore_funcs:
                    continue
                recipes.append(
                    Recipe(prefix,
                           func_name,
                           bases,
                           func.get("implements", ""),
                           func.get("relates", "")))
        return HeaderFile(self._path, [], recipes)

    @classmethod
    def parse_h(cls, ignore_files, ignore_funcs) -> List[HeaderFile]:
        """Parse existing header file."""
        files = [_ for _ in _clib_path.glob("*.h")
                 if _.name not in ignore_files + _clib_ignore]
        return [cls(_, ignore_funcs.get(_.name, []))._parse_h() for _ in files]

    def _parse_h(self) -> HeaderFile:
        ct = self._path.read_text()

        def parse_with_doxygen(text):
            # a primitive doxygen parser for existing Cantera CLib header files
            regex = re.compile((
                r"(?P<blank>\n\s*\n)|"  # blank line
                r"(?P<head>(?=//! )[^\n]*)|"  # leading doxygen comment
                r"(?P<func>(?=CANTERA_CAPI)[^;]*;)|"  # CLib function
                r"(?P<tail>(?=//!< )[^\n]*)"))  # trailing doxygen comment
            matches = re.finditer(regex, text)

            pairs = []
            function = None
            comments = []
            for m in matches:
                if function:
                    if m.group("func"):
                        # new function: flush buffers
                        pairs.append((function, "\n".join(comments)))
                        comments = []
                        function = m.group("func")
                    elif m.group("head"):
                        # new heading comment: flush buffers
                        pairs.append((function, "\n".join(comments)))
                        comments = [m.group("head")]
                        function = None
                    elif m.group("tail"):
                        # trailing comment: append to buffer
                        comments.append(m.group("tail"))
                    else:  # m.group("blank"):
                        # blank line: flush buffers
                        pairs.append((function, "\n".join(comments)))
                        comments = []
                        function = None
                else:
                    if m.group("blank"):
                        # blank line: clear buffer
                        comments = []
                    elif m.group("head"):
                        # new heading comment: buffer
                        comments.append(m.group("head"))
                    elif m.group("func"):
                        # new function: buffer
                        function = m.group("func")
            return pairs

        c_functions = parse_with_doxygen(ct)

        if not c_functions:
            return

        parsed = map(self._parse_func, c_functions)

        logger.info(f"  parsing {self._path.name}")
        if self._ignore_funcs:
            logger.info(f"    ignoring {self._ignore_funcs}")

        parsed = [f for f in parsed if f.name not in self._ignore_funcs]

        if not parsed:
            return

        return HeaderFile(self._path, parsed)


def doxygen_func(tag: str, text: str) -> str:
    """Extract function signature from doxygen tag if it exists."""
    regex = re.compile(rf"(?<={tag} ).*[^\(\n]|((?<={tag} )(.*?)\))")
    matched = list(re.finditer(regex, text))
    if not matched:
        return ""
    if len(matched) > 1:
        msg = f"Found more than one {tag!r} annotation; returning first."
        logging.warning(msg)
        signatures = '\n  - '.join([""] + [_[0] for _ in matched])
        logging.debug("Found instances:%s", signatures)
    return matched[0][0]
