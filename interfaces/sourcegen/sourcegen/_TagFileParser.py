# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from pathlib import Path
import re
import logging
from typing import Dict, List, Union

from ._dataclasses import Func, AnnotatedFunc, ArgList


logger = logging.getLogger(__name__)

_tag_path = Path(__file__).parent.joinpath("../../../build/doc/").resolve()

class TagFileParser:
    """Class handling contents of doxygen tag file."""

    def _parse_doxyfile(self, class_crosswalk: List[str]):
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
            msg = f"Missing {kind!r} compound(s):\n    {missing!r}\nusing regex "
            msg += f"{regex}. Continuing with remaining compounds: \n    {found!r}"
            logger.error(msg)

        # Parse content of namespace Cantera
        namespace = xml_compounds("namespace", ["Cantera"])["Cantera"]
        qualified_names = xml_tags("class", namespace, suffix='kind="class"')
        class_names = [_.split(":")[-1] for _ in qualified_names]

        # Handle exceptions for unknown/undocumented classes
        unknown = set(class_crosswalk) - set(class_names)
        if unknown:
            unknown = "', '".join(unknown)
            logger.warning(
                "Class(es) in configuration file are missing from tag file: "
                f"'{unknown}'")

        # Parse content of classes that are specified by the configuration file
        class_names = set(class_crosswalk) & set(class_names)
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
            msg = (f"Tag file does not exist at expected location:\n    {tag_file}\n"
                   "Run 'scons doxygen' to generate.")
            logger.critical(msg)
            sys.exit(1)

        with tag_file.open() as fid:
            self._doxygen_tags = fid.read()

        self._parse_doxyfile(class_crosswalk)

    def annotated_func(self, parsed: Func) -> Union[AnnotatedFunc, None]:
        """Match function with doxygen tag information."""
        comments = parsed.annotations

        def doxygen_func(tag: str, text: str) -> Union[str, None]:
            regex = re.compile(r"(?<={0} ).*[^\(\n]|((?<={0} )(.*?)\))".format(tag))
            matched = list(re.finditer(regex, text))
            if not matched:
                return None
            if len(matched) > 1:
                msg = f"Found more than one {tag!r} annotation; returning first."
                logging.warning(msg)
            return matched[0][0]

        implements = doxygen_func("@implements", comments)
        if not implements:
            return None

        cxx_func = implements.split("(")[0]
        if cxx_func not in self._known:
            logger.error(f"Did not find {cxx_func!r} in tag file.")
            return None
        ix = 0
        if len(self._known[cxx_func]) > 1:
            # Disambiguate functions with same name
            args = re.findall(re.compile(r'(?<=\().*(?=\))'), implements)
            if not args:
                logger.error(
                    f"Need argument list to disambiguate {implements!r}: skipping.")
                return None
            args = f"({args[0]})"
            ix = -1
            for i, xml in enumerate(self._known[cxx_func]):
                arglist = ArgList.from_xml(xml_tags("arglist", xml)[0])
                if args == arglist.short_str():
                    ix = i
                    break
            if ix < 0:
                logger.error(
                    f"Unable to match {cxx_func!r} to known functions: skipping.")
                return None
        xml = self._known[cxx_func][ix]

        return AnnotatedFunc(*parsed,
                             implements,
                             doxygen_func("@relates", comments),
                             xml_tags("type", xml)[0],
                             xml_tags("name", xml)[0],
                             xml_tags("arglist", xml)[0],
                             xml_tags("anchorfile", xml)[0].replace(".html", ".xml"),
                             xml_tags("anchor", xml)[0])


def xml_tags(tag: str, text: str, suffix: str="") -> Union[str, None]:
    if suffix:
        suffix = f" {suffix.strip()}"
    regex = re.compile(r'(?<=<{0}{1}>)(.*?)(?=</{0}>)'.format(tag, suffix))
    matched = re.findall(regex, text)
    if not matched:
        blanks = text.split("\n")[-1].split("<")[0]
        msg = f"Could not extract {tag!r} from:\n{blanks}{text}\n"
        msg += f"using regex: {regex}"
        logger.error(msg)
        return None
    return matched
