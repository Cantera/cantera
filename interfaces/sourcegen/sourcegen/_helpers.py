# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import fields
import inspect
from pathlib import Path
import textwrap
try:
    from ruamel import yaml
except ImportError:
    import ruamel_yaml as yaml


def read_config(config_file: Path) -> dict:
    """Read YAML configuration file."""
    if config_file.is_file():
        with config_file.open("r", encoding="utf-8") as config_file:
            reader = yaml.YAML(typ="safe")
            return reader.load(config_file)
    else:
        return {}


def with_unpack_iter(cls: type) -> type:
    cls.__iter__ = lambda self: (getattr(self, f.name) for f in fields(self))

    return cls


def hanging_text(text: str, spaces: int) -> str:
    ret = ("\n" + " "*spaces).join(text.split("\n"))
    return "\n".join([line.rstrip() for line in ret.split("\n")])


def normalize_indent(code: str) -> str:
    code = textwrap.dedent(code).strip()

    call_line = inspect.stack()[1].code_context[0]

    indent = len(call_line) - len(call_line.lstrip())

    # If called inside a string interpolation, indent to the rest of the block.
    # Look for the opening brace for the interpolation, which isn"t perfect, but works.
    # This will fire for lines such as “        {normalize(my_str)}”
    if call_line[indent] == "{":
        code = textwrap.indent(code, " " * indent)

        code = code[indent:]

    return code


def get_preamble() -> str:
    return Path(__file__).parent.joinpath("preamble.txt").read_text("utf-8").strip()
