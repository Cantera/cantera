# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import re
import warnings
from typing import Dict, List, Union

from ._dataclasses import Func, AnnotatedFunc


_tag_path = Path(__file__).parent.joinpath("../../../build/doc/").resolve()

class TagFileParser:
    """Class handling contents of doxygen tag file."""

    def _parse_doxyfile(self, class_crosswalk: Dict[str, str]):
        """Retrieve class and function information from Cantera namespace."""

        def xml_compounds(kind: str, name: str="") -> List[str]:
            regex = re.compile(r'<compound kind="{0}"[\s\S]*?</compound>'.format(kind))
            if not name:
                return re.findall(regex, self._doxygen_tags)
            for compound in re.findall(regex, self._doxygen_tags):
                compound_name = xml_tags("name", compound)[0]
                if compound_name == name:
                    return [compound]
            msg = f"Tag file does not contain compound {kind!r} with name {name!r}."
            raise ValueError(msg)

        namespace = xml_compounds("namespace", "Cantera")[0]

        # Get class names and handle exceptions for unknown/undocumented classes
        qualified_names = xml_tags("class", namespace, suffix='kind="class"')
        class_names = [_.split(":")[-1] for _ in qualified_names]
        unknown = set(class_crosswalk.values()) - set(class_names)
        if unknown:
            unknown = '", "'.join(unknown)
            msg = f"Unknown/undocumented class(es) in configuration file: {unknown!r}"
            warnings.warn(msg)

        # Parse content of classes that are specified by the configuration file
        class_names = set(class_crosswalk.values()) & set(class_names)
        regex = re.compile(r'<compound kind="class">[\s\S]*?</compound>')
        classes = {}
        regex_name = re.compile(r'(?<=<name>)(.*?)(?=</name>)')
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

        self._classes = classes
        self._functions = functions

    def __init__(self, class_crosswalk: Dict[str, str]) -> None:
        tag_file = _tag_path / "Cantera.tag"
        if not tag_file.exists():
            raise FileNotFoundError(
                f"Tag file does not exist at expected location:\n    {tag_file}\n"
                "Run 'scons doxygen' to generate.")

        with tag_file.open() as fid:
            self._doxygen_tags = fid.read()

        self._parse_doxyfile(class_crosswalk)

    @staticmethod
    def from_parsed(parsed_config_file: dict):
        return TagFileParser(parsed_config_file["class_crosswalk"])

    def annotated_func(self, parsed: Func) -> Union[AnnotatedFunc, None]:
        """Match function with doxygen tag information."""
        ret_type, name, params, comments = parsed

        def doxygen_func(tag: str, text: str) -> Union[str, None]:
            regex = re.compile(r"(?<={0} ).*[^\(\n]|((?<={0} )(.*?)\))".format(tag))
            matched = list(re.finditer(regex, text))
            if not matched:
                return None
            if len(matched) > 1:
                msg = f"Found more than one {tag!r} annotation."
                raise RuntimeError(msg)
            return matched[0][0]

        implements = doxygen_func("@implements", comments)
        if not implements:
            return None
        relates = doxygen_func("@relates", comments)

        cxx_func = implements.split("(")[0]
        if "::" in cxx_func:
            # A class
            return AnnotatedFunc(ret_type, name, params, comments, implements, relates, "", "", "", "", "")

        if cxx_func not in self._functions:
            msg = f"Did not find {cxx_func!r} in tag file."
            warnings.warn(msg)
            return None
        if len(self._functions[cxx_func]) > 1:
            msg = f"Ambiguous function {cxx_func!r}"
            warnings.warn(msg)
            return AnnotatedFunc(ret_type, name, params, comments, implements, relates, "", "", "", "", "")
        xml = self._functions[cxx_func][0]

        # A function
        return AnnotatedFunc(*parsed,
                             implements,
                             relates,
                             xml_tags("type", xml)[0],
                             xml_tags("name", xml)[0],
                             xml_tags("anchorfile", xml)[0].replace(".html", ".xml"),
                             xml_tags("anchor", xml)[0],
                             xml_tags("arglist", xml)[0])

def xml_tags(tag: str, text: str, suffix: str="") -> Union[str, None]:
    if suffix:
        suffix = f" {suffix.strip()}"
    regex = re.compile(r'(?<=<{0}{1}>)(.*?)(?=</{0}>)'.format(tag, suffix))
    matched = re.findall(regex, text)
    if not matched:
        blanks = text.split("\n")[-1].split("<")[0]
        msg = f"Could not extract {tag!r} from:\n{blanks}{text}\n"
        msg += f"using regex: {regex}"
        raise RuntimeError(msg)
    return matched
