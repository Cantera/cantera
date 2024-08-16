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

        def xml_compounds(kind: str, names: List[str]) -> Dict[str,str]:
            regex = re.compile(r'<compound kind="{0}"[\s\S]*?</compound>'.format(kind))
            found = []
            compounds = {}
            for compound in re.findall(regex, self._doxygen_tags):
                qualified_name = xml_tags("name", compound)[0]
                compound_name = qualified_name.split(":")[-1]
                if compound_name in names:
                    found.append(compound_name)
                    compounds[compound_name] = compound
                    if not set(names) - set(found):
                        return compounds
            missing = '", "'.join(set(names) - set(found))
            msg = f"Missing {kind!r} compound(s):\n    {missing!r}\n"
            msg += f"using regex: {regex}"
            raise ValueError(msg)

        # Parse content of namespace Cantera
        namespace = xml_compounds("namespace", ["Cantera"])["Cantera"]
        qualified_names = xml_tags("class", namespace, suffix='kind="class"')
        class_names = [_.split(":")[-1] for _ in qualified_names]

        # Handle exceptions for unknown/undocumented classes
        unknown = set(class_crosswalk.values()) - set(class_names)
        if unknown:
            unknown = '", "'.join(unknown)
            msg = f"Unknown/undocumented class(es) in configuration file: {unknown!r}"
            warnings.warn(msg)

        # Parse content of classes that are specified by the configuration file
        class_names = set(class_crosswalk.values()) & set(class_names)
        classes = xml_compounds("class", class_names)

        def xml_members(kind: str, text: str, prefix="") -> Dict[str, str]:
            regex = re.compile(r'<member kind="{0}"[\s\S]*?</member>'.format(kind))
            functions = {}
            for func in re.findall(regex, text):
                func_name = f'{prefix}{xml_tags("name", func)[0]}'
                if func_name in functions:
                    functions[func_name].append(func)
                else:
                    functions[func_name] = [func]
            return functions

        # Get known functions from namespace and methods from classes
        self._known = xml_members("function", namespace)
        for name, cls in classes.items():
            prefix = f"{name}::"
            self._known.update(xml_members("function", cls, prefix))

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

        def shorten_arglist(arglist: str) -> str:
            # Only keep contents between parentheses and replace HTML entities
            ret = re.findall(re.compile(r'(?<=\().*(?=\))'), arglist)[0]
            replacements = [(" &amp;", "& "), ("&lt; ", "<"), (" &gt;", ">")]
            for rep in replacements:
                ret = ret.replace(*rep)
            # Remove parameters names and default values from argument string
            return ", ".join([" ".join(_.split()[:-1]) for _ in ret.split(",")])

        cxx_func = implements.split("(")[0]
        if cxx_func not in self._known:
            msg = f"Did not find {cxx_func!r} in tag file."
            warnings.warn(msg)
            return None
        ix = 0
        if len(self._known[cxx_func]) > 1:
            # Disambiguate functions with same name
            args = re.findall(re.compile(r'(?<=\().*(?=\))'), implements)
            if not args:
                msg = f"Need argument list to disambiguate {implements!r}"
                raise RuntimeError(msg)
            args = args[0]
            ix = -1
            for i, xml in enumerate(self._known[cxx_func]):
                arglist = xml_tags("arglist", xml)[0]
                if args == shorten_arglist(arglist):
                    ix = i
                    break
            if ix < 0:
                msg = f"Unable to match {cxx_func!r} to known function."
                raise RuntimeError(msg)
        xml = self._known[cxx_func][ix]

        return AnnotatedFunc(*parsed,
                             implements,
                             doxygen_func("@relates", comments),
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
