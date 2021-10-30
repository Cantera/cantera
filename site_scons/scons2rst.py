#!/usr/bin/env python3
# encoding: utf-8

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
scons2rst.py: Extract SCons configuration options from SConstruct
"""

import ast
import re
import argparse
import pathlib


def deblank(string):
    """
    Remove whitespace before and after line breaks
    """
    out = [s.strip() for s in string.split("\n")]
    if not len(out[-1]):
        out = out[:-1]
    return "\n".join(out)


def convert(item):
    """
    Convert AST objects to human-readable strings (or list of strings)
    """
    if isinstance(item, list):
        return [convert(sub) for sub in item]

    if isinstance(item, (ast.Tuple, ast.List)):
        return [convert(sub) for sub in item.elts]

    if isinstance(item, ast.Str):
        return deblank(item.s)

    if isinstance(item, (ast.NameConstant)):
        if isinstance(item.value, bool):
            return "True" if item.value else "False"
        raise ValueError(f"Unable to process constant '{item}'")

    if isinstance(item, ast.Attribute):
        if isinstance(item.value, ast.Name):
            return f"{item.value.id}.{item.attr}"
        raise ValueError(f"Unable to process attribute '{item}'")

    if isinstance(item, ast.Subscript):
        if isinstance(item.value, ast.Name) and isinstance(item.slice.value, ast.Str):
                return f"{item.value.id}[{item.slice.value.s}]"
        raise ValueError(f"Unable to process subscript '{item}'")

    if isinstance(item, ast.Name):
        return item.id

    raise ValueError(f"Unable to process item of type '{type(item)}'")


class Option:
    """
    Object corresponding to SCons configuration option
    """

    def __init__(self, item):
        if isinstance(item, (ast.Tuple, ast.List)):
            self.func = None
            self._assign(convert(item))
            return

        if isinstance(item, ast.Call):
            self.func = item.func.id
            self._assign(convert(item.args))
            return

        raise ValueError(f"Unable to process item of type '{type(item)}'")

    def _assign(self, args):
        self.name = args[0]
        self.description = args[1]
        self.default = args[2]
        self.choices = args[3] if len(args) > 3 else None

    @property
    def args(self):
        """
        Arguments of the SCons configuration option
        """
        if self.choices is None:
            return self.name, self.description, self.default
        return self.name, self.description, self.default, self.choices

    def to_rest(self, dev=False, indent=3):
        """
        Convert option to restructured text
        """
        tag = self.name.replace("_", "-").lower()
        if dev:
            tag += "-dev"

        if self.func == "PathVariable":
            choices = f"``path/to/{self.name}``"
            default = self.default
        elif isinstance(self.choices, list):
            choices = self.choices
            for yes_no in ["n", "y"]:
                if yes_no in choices:
                    # ensure correct order
                    choices.remove(yes_no)
                    choices = [yes_no] + choices
            choices = " | ".join([f"``{c}``" for c in choices])
            default = self.default
        elif self.default in ["True", "False"]:
            choices = "``yes`` | ``no``"
            default = "yes" if self.default == "True" else "no"
        elif self.func == "BoolVariable":
            choices = "``yes`` | ``no``"
            default = self.default
        else:
            choices = "``string``"
            default = self.default

        if default.startswith("defaults.") or default.startswith("sys."):
            default = "''"
        else:
            default = f"'{default}'"

        # assemble description
        tab = " " * indent
        linebreak = "\n" + tab
        description = linebreak.join(self.description.split("\n"))
        pat = r'"([a-zA-Z0-9\-\+_.,: =/\'\\]+)"'
        double_quoted = []
        for item in re.findall(pat, description):
            # enclose double-quoted strings in '``'
            found = f'"{item}"'
            double_quoted += [found]
            replacement = f"``{found}``"
            description = description.replace(found, replacement)
        pat = r"\'([a-zA-Z0-9\-\+_.,:=/\\]+)\'"
        for item in re.findall(pat, description):
            # replace "'" for single-quoted words by '``'; do not replace "'" when
            # whitespace is enclosed or if word is part of double-quoted string
            if any([item in dq for dq in double_quoted]):
                continue
            found = f"'{item}'"
            replacement = found.replace("'", "``")
            description = description.replace(found, replacement)
        pat = r"\*([a-zA-Z0-9\-\+_., :=/\'\\]+)"
        asterisks = re.findall(pat, description)
        if len(asterisks) == 1:
            # catch unbalanced '*', for example in '*nix'
            found = f"*{asterisks[0]}"
            replacement = f"\{found}"
            description = description.replace(found, replacement)

        # assemble output
        out = f".. _{tag}:\n\n"
        out += f"*  ``{self.name}``: [ {choices} ]\n"
        out += f"{tab}{description}\n\n"
        out += f"{tab}- default: ``{default}``\n"

        return out

    def __repr__(self):
        if self.func is None:
            return f"{self.args}"
        return f"{self.func}{self.args}"


def parse(input_file, output_file, dev):
    """Parse SConstruct and extract configuration"""

    with open(input_file, "r") as fid:
        code = fid.read()

    tree = ast.parse(code)

    names = ["windows_compiler_options", "compiler_options", "config_options"]
    options = {name: [] for name in names}

    for node in ast.walk(tree):

        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name) \
                and len(node.targets) == 1:
            # option is assigned as list
            name = node.targets[0].id
            if name not in names:
                continue

            if isinstance(node.value, ast.List):
                for item in node.value.elts:
                    options[name].append(Option(item))
                continue

            if isinstance(node.value, ast.Call):
                continue

            raise ValueError(f"Unable to process node {node.value} (name='{name}')")

        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and isinstance(node.func.value, ast.Name):
            # option is appended to list via 'append' or 'extend'
            if node.func.attr not in ["append", "extend"]:
                continue

            name = node.func.value.id
            if name not in names:
                continue

            for arg in node.args:
                if isinstance(arg, ast.List):
                    for item in arg.elts:
                        options[name].append(Option(item))
                    continue

                if isinstance(arg, ast.Call):
                    options[name].append(Option(arg))
                    continue

                raise ValueError(f"Unable to process node {node}")

    with open(output_file, "w+") as fid:
        for config in names:
            for option in options[config]:
                fid.write(option.to_rest(dev=dev))
                fid.write("\n")

    print(f"Done writing output to '{output_file}'.")

def main():
    parser = argparse.ArgumentParser(
        description="Extract configuration options from SConstruct",
    )
    parser.add_argument(
        "input",
        nargs="?",
        help="The input SConstruct file. Optional.",
        default="SConstruct")
    parser.add_argument(
        "output",
        nargs="?",
        help="The output ReST filename. Optional.",
        default="config-options")
    parser.add_argument(
        "--dev",
        action="store_true",
        default=False,
        help="Append '-dev' to filename and Sphinx references.")

    args = parser.parse_args()
    input_file = pathlib.Path(args.input)
    output_file = args.output
    if args.dev:
        output_file += "-dev"
    output_file = pathlib.Path(output_file).with_suffix(".rst")

    parse(input_file, output_file, args.dev)


if __name__ == "__main__":
    main()
