# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging
import re
from typing import List

from ._dataclasses import HeaderFile, Func, Param


logger = logging.getLogger()

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

    def parse(self) -> HeaderFile:
        ct = self._path.read_text()

        def parse_with_doxygen(text):
            # a primitive doxygen parser for Cantera CLib header files
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
