# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from dataclasses import fields
import inspect
from pathlib import Path
import textwrap
import re
try:
    from ruamel import yaml
except ImportError:
    import ruamel_yaml as yaml


def read_config(config_file: Path) -> dict:
    if config_file.is_file():
        with config_file.open("r", encoding="utf-8") as config_file:
            reader = yaml.YAML(typ="safe")
            return reader.load(config_file)
    else:
        return {}


def with_unpack_iter(cls: type) -> type:
    cls.__iter__ = lambda self: (getattr(self, f.name) for f in fields(self))

    return cls


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
    return Path(__file__).parent.joinpath("preamble.txt").read_text("utf-8")


def xml_tag(tag: str, text: str, suffix: str="", index=0) -> str:
    """Extract content enclosed between XML tags, optionally skipping matches."""
    if suffix:
        suffix = f" {suffix.strip()}"
    regex = re.compile(rf'(?<=<{tag}{suffix}>)(.*?)(?=</{tag}>)')
    match = re.findall(regex, text)
    if index >= len(match):
        return ""  # not enough matches found
    return match[index]


def split_arglist(arglist: str) -> tuple:
    """Split C++ argument list into text within parentheses and suffix."""
    arglist = arglist.strip()
    suffix = arglist[arglist.rfind(")") + 1:]
    arglist = re.findall(re.compile(r'(?<=\().*(?=\))'), arglist)[0]
    return arglist, suffix.strip()
