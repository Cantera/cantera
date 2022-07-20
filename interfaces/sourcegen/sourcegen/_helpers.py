# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import os
import textwrap

def normalize(code: str, indent: int = 0, trim_first: bool = False):
    code = textwrap.dedent(code).strip()

    if indent:
        code = textwrap.indent(code, ' ' * indent)

        if trim_first:
            code = code[indent:]

    return code

def get_preamble():
    with open(os.path.dirname(__file__) + '/preamble.txt', 'r') as preamble_file:
        return preamble_file.read()
