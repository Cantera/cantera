# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
from typing import List, Dict, Union
import warnings
import re

from ._Config import Config

from .._dataclasses import Func, HeaderFile, AnnotatedFunc
from .._SourceGenerator import SourceGenerator


_tag_path = Path(__file__).parent.joinpath("../../../../build/doc/").resolve()

class DoxygenSourceGenerator(SourceGenerator):
    """The SourceGenerator for referencing CLib functions to doxygen information."""

    def _parse_doxyfile(self):
        # Get contents of Cantera namespace
        regex = re.compile(r'<compound kind="namespace"[\s\S]*?</compound>')
        regex_name = re.compile(r'(?<=<name>)(.*?)(?=</name>)')
        namespace = None
        for ns in re.findall(regex, self._doxygen_tags):
            names = re.findall(regex_name, ns)
            if names and "Cantera" == names[0]:
                namespace = ns
                break
        if not namespace:
            raise ValueError("Invalid tag file does not contain namespace 'Cantera'.")

        # Get class names and handle exceptions for unknown/undocumented classes
        regex = re.compile(r'(?<=<class kind="class">Cantera::)(.*?)(?=</class>)')
        class_names = re.findall(regex, namespace)
        unknown = set(self._config.class_crosswalk.values()) - set(class_names)
        if unknown:
            unknown = '", "'.join(unknown)
            msg = f"Unknown/undocumented class(es) in configuration file: {unknown!r}"
            warnings.warn(msg)

        # Parse content of classes that are specified by the configuration file
        class_names = set(self._config.class_crosswalk.values()) & set(class_names)
        qualified_names = [f"Cantera::{_}" for _ in class_names]
        regex = re.compile(r'<compound kind="class">[\s\S]*?</compound>')
        classes = {}
        for cl in re.findall(regex, self._doxygen_tags):
            names = re.findall(regex_name, cl)
            if names and names[0] in qualified_names:
                classes[names[0].split(":")[-1]] = cl

        # Get functions defined in Cantera namespace
        regex = re.compile(r'<member kind="function">[\s\S]*?</member>')
        functions = {}
        for fn in re.findall(regex, namespace):
            names = re.findall(regex_name, fn)
            if names:
                name = names[0]
                if name in functions:
                    functions[name].append(fn)
                else:
                    functions[name] = [fn]

        return classes, functions

    def _convert_func(self, parsed: Func) -> Union[AnnotatedFunc, None]:
        ret_type, name, params, comments = parsed
        if not comments or "@implements" not in comments:
            return None
        regex = re.compile(r"(?<=@implements ).*[^\(\n]|((?<=@implements )(.*?)\))")
        found = list(re.finditer(regex, comments))
        if not found:
            raise RuntimeError("This needs debugging.")

        implemented = found[0][0]
        return AnnotatedFunc(ret_type, name, params, comments, implemented, "", "", "", "", "")

    def _find_anchor(self, func):
        if "::" in func.implements:
            # A class
            print(f"{func.name:<18}: {func.implements}")
        elif func.implements not in self._functions:
            msg = f"Unable to find {func.implements!r} in doxygen tag file"
            warnings.warn(msg)
            return
        else:
            # A function
            if len(self._functions[func.implements]) > 1:
                raise RuntimeError(f"Ambiguous function {func.implemens!r}")
            xml = self._functions[func.implements][0]
            regex = re.compile(r'(?<=<anchorfile>)(.*?)(?=</anchorfile>)')
            anchor_file = re.findall(regex, xml)
            regex = re.compile(r'(?<=<anchor>)(.*?)(?=</anchor>)')
            anchor = re.findall(regex, xml)
            if anchor:
                print(f"{func.name:<18}: {func.implements:<20} "
                      f"{anchor_file[0].replace('.html', '.xml')} ({anchor[0]})")

    def __init__(self, out_dir: str, config: dict):
        if out_dir:
            raise ValueError("Doxygen generator does not use output folder.")
        tag_file = _tag_path / "Cantera.tag"
        if not tag_file.exists():
            raise FileNotFoundError(
                f"Tag file does not exist at expected location:\n    {tag_file}\n"
                "Run 'scons doxygen' to generate.")

        with tag_file.open() as fid:
            self._doxygen_tags = fid.read()

        # use the typed config
        self._config = Config.from_parsed(config)
        self._classes, self._functions = self._parse_doxyfile()

    def generate_source(self, headers_files: List[HeaderFile]):
        documented_funcs: Dict[str, AnnotatedFunc] = {}
        undocumented_funcs: List[str] = []

        for header_file in headers_files:
            doxy_funcs = list(map(self._convert_func, header_file.funcs))
            documented_funcs.update((f.name, f) for f in doxy_funcs if f is not None)
            undocumented_funcs.extend([_ for _ in doxy_funcs if _ is None])

        for func in documented_funcs.values():
            self._find_anchor(func)
            # print(f"{func.name}: {func.implements}")
