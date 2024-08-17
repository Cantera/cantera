# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from pathlib import Path
import re
import logging
from typing import Dict, List, Union

from ._dataclasses import ArgList, TagInfo
from ._helpers import xml_tag


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
                qualified_name = xml_tag("name", compound)
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

        # Parse content of namespace Cantera
        namespace = xml_compounds("namespace", ["Cantera"])["Cantera"]
        qualified_names = xml_tags("class", namespace, suffix='kind="class"')
        class_names = [_.split(":")[-1] for _ in qualified_names]

        # Handle exceptions for unknown/undocumented classes
        unknown = set(class_crosswalk) - set(class_names)
        if unknown:
            unknown = "', '".join(unknown)
            logger.warning(
                "Class(es) in configuration file are missing from tag file: '%s'",
                unknown)

        # Parse content of classes that are specified by the configuration file
        class_names = set(class_crosswalk) & set(class_names)
        classes = xml_compounds("class", class_names)

        def xml_members(kind: str, text: str, prefix="") -> Dict[str, str]:
            regex = re.compile(rf'<member kind="{kind}"[\s\S]*?</member>')
            functions = {}
            for func in re.findall(regex, text):
                func_name = f'{prefix}{xml_tag("name", func)}'
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

    def tag_info(self, func_string: str) -> TagInfo:
        """Look up tag information based on (partial) function signature."""
        cxx_func = func_string.split("(")[0].split(" ")[-1]
        if cxx_func not in self._known:
            logger.error("Did not find '%s' in tag file.", cxx_func)
            return None
        ix = 0
        if len(self._known[cxx_func]) > 1:
            # Disambiguate functions with same name
            args = re.findall(re.compile(r'(?<=\().*(?=\))'), func_string)
            if not args:
                logger.error("Need argument list to disambiguate '%s': skipping.",
                             func_string)
                return TagInfo()
            args = f"({args[0]})"
            ix = -1
            for i, xml in enumerate(self._known[cxx_func]):
                arglist = ArgList.from_xml(xml_tag("arglist", xml))
                if args == arglist.short_str():
                    ix = i
                    break
            if ix < 0:
                logger.error("Unable to match '%s' to known functions: skipping.",
                             func_string)
                return TagInfo()
        return TagInfo.from_xml(self._known[cxx_func][ix])
