# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import dataclass
from pathlib import Path
import re
from ._helpers import with_unpack_iter


@dataclass(frozen=True)
@with_unpack_iter
class Param:
    """ Represents a function parameter """

    p_type: str
    name: str


    @staticmethod
    def parse(c_param: str):
        parts = c_param.strip().rsplit(' ', 1)
        if len(parts) == 2:
            return Param(*parts)


@dataclass(frozen=True)
@with_unpack_iter
class Func:
    """ Represents a function parsed from a C header file """

    ret_type: str
    name: str
    params: list[Param]


    @staticmethod
    def parse(c_func: str):
        lparen = c_func.index('(')
        rparen = c_func.index(')')
        front = c_func[0:lparen].split()

        params = list(map(Param.parse,
            filter(None, c_func[lparen+1:rparen].split(','))))

        ret_type = front[-2]
        name = front[-1]
        return Func(ret_type, name, params)


@dataclass(frozen=True)
@with_unpack_iter
class HeaderFile:
    """ Represents information about a parsed C header file """

    path: Path
    funcs: list[Func]


    @staticmethod
    def parse(file: Path, ignore: list[str]):
        ct = file.read_text()

        matches = re.finditer(r'CANTERA_CAPI.*?;', ct, re.DOTALL)
        c_functions = [re.sub(r'\s+', ' ', m.group()) for m in matches]

        if not c_functions:
            return

        parsed = map(Func.parse, c_functions)

        if ignore:
            print(f'  ignoring ' + str(ignore))

        parsed = [f for f in parsed if f.name not in ignore]

        if not parsed:
            return

        return HeaderFile(file, parsed)
