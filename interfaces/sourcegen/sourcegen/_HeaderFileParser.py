# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import re
from typing import List

from ._dataclasses import HeaderFile, Func, Param


class HeaderFileParser:

    @staticmethod
    def _parse_param(c_param: str) -> Param:
        parts = c_param.strip().rsplit(" ", 1)
        if len(parts) == 2:
            return Param(*parts)

    @classmethod
    def _parse_func(cls, c_func: str) -> Func:
        lparen = c_func.index("(")
        rparen = c_func.index(")")
        front = c_func[0:lparen].split()

        params = list(map(cls._parse_param,
            filter(None, c_func[(lparen + 1):rparen].split(","))))

        ret_type = front[-2]
        name = front[-1]
        return Func(ret_type, name, params)

    def __init__(self, path: Path, ignore_funcs: List[str] = None):
        self._path = path
        self._ignore_funcs = ignore_funcs

    def parse(self) -> HeaderFile:
        ct = self._path.read_text()

        matches = re.finditer(r"CANTERA_CAPI.*?;", ct, re.DOTALL)
        c_functions = [re.sub(r"\s+", " ", m.group()) for m in matches]

        if not c_functions:
            return

        parsed = map(self._parse_func, c_functions)

        print(f"  parsing " + self._path.name)
        if self._ignore_funcs:
            print(f"    ignoring " + str(self._ignore_funcs))

        parsed = [f for f in parsed if f.name not in self._ignore_funcs]

        if not parsed:
            return

        return HeaderFile(self._path, parsed)
