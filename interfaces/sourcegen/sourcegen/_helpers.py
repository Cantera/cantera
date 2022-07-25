# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import inspect
import os
import textwrap

def normalize_indent(code: str):
    code = textwrap.dedent(code).strip()

    indent = 0

    call_line = inspect.stack()[1].code_context[0]

    while call_line[indent] == ' ':
        indent += 1

    # If called inside a string interpolation, indent to the rest of the block.
    # Look for the opening brace for the interpolation, which isn't perfect, but works.
    # This will fire for lines such as “        {normalize(my_str)}”
    if call_line[indent] == '{':
        code = textwrap.indent(code, ' ' * indent)

        code = code[indent:]

    return code

def get_preamble():
    with open(os.path.dirname(__file__) + '/preamble.txt', 'r') as preamble_file:
        return preamble_file.read()
