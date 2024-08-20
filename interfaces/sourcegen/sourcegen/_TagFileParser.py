"""Parser for tag files and XML generated by doxygen."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from pathlib import Path
import re
import logging
from dataclasses import dataclass
from typing import Dict, List

from ._dataclasses import ArgList, Param
from ._helpers import with_unpack_iter


logger = logging.getLogger(__name__)

_tag_path = Path(__file__).parent.joinpath("../../../build/doc/").resolve()
_xml_path = _tag_path / "doxygen" / "xml"


@dataclass(frozen=True)
@with_unpack_iter
class TagInfo:
    """Represents information parsed from a doxygen tag file."""

    base: str = ""  #: qualified scope (skipping Cantera namespace)
    type: str = ""  #: return type
    name: str = ""  #: function name
    arglist: str = ""  #: function argument list
    anchorfile: str = ""  #: doxygen anchor file
    anchor: str = ""  #: doxygen anchor

    @staticmethod
    def from_xml(qualified_name, xml):
        """Create tag information based on XML data."""
        base = ""
        if "::" in qualified_name:
            base = qualified_name.split("::", 1)[0]
        return TagInfo(base,
                       xml_tag("type", xml),
                       xml_tag("name", xml),
                       xml_tag("arglist", xml),
                       xml_tag("anchorfile", xml).replace(".html", ".xml"),
                       xml_tag("anchor", xml))

    def __bool__(self):
        return all([self.type, self.name, self.arglist, self.anchorfile, self.anchor])

    @property
    def signature(self):
        """Generate function signature based on tag information."""
        ret = f"{self.type} {self.name}{self.arglist}"
        replacements = [(" &amp;", "& "), ("&lt; ", "<"), (" &gt;", ">")]
        for rep in replacements:
            ret = ret.replace(*rep)
        return ret

    @property
    def id(self):
        """Generate doxygen id."""
        return f"{self.anchorfile.replace('.xml', '')}_1{self.anchor}"

    @property
    def qualified_name(self):
        """Return qualified name."""
        if self.base:
            return f"{self.base}::{self.name}"
        return self.name


@dataclass(frozen=True)
@with_unpack_iter
class TagDetails(TagInfo):
    """Create tag information based on XML data."""

    briefdescription: str = ""
    location: str = ""
    parameterlist: List[Param] = None


class TagFileParser:
    """Class handling contents of doxygen tag file."""

    def _parse_doxyfile(self, bases: List[str]):
        """Retrieve class and function information from Cantera namespace."""

        def xml_compounds(kind: str, names: List[str]) -> Dict[str,str]:
            regex = re.compile(rf'<compound kind="{kind}"[\s\S]*?</compound>')
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

        # Parse content of namespace Cantera
        namespace = xml_compounds("namespace", ["Cantera"])["Cantera"]
        qualified_names = xml_tags("class", namespace, suffix='kind="class"')
        class_names = [_.split(":")[-1] for _ in qualified_names]

        # Handle exceptions for unknown/undocumented classes
        unknown = set(bases) - set(class_names)
        if unknown:
            unknown = "', '".join(unknown)
            logger.warning(
                "Class(es) in configuration file are missing from tag file: '%s'",
                unknown)

        # Parse content of classes that are specified by the configuration file
        class_names = set(bases) & set(class_names)
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

    def __init__(self, bases: Dict[str, str]) -> None:
        tag_file = _tag_path / "Cantera.tag"
        if not tag_file.exists():
            msg = (f"Tag file does not exist at expected location:\n    {tag_file}\n"
                   "Run 'scons doxygen' to generate.")
            logger.critical(msg)
            sys.exit(1)

        with tag_file.open() as fid:
            self._doxygen_tags = fid.read()

        logging.info("Parsing doxygen tags...")
        self._parse_doxyfile(bases)

    def tag_info(self, func_string: str) -> TagInfo:
        """Look up tag information based on (partial) function signature."""
        cxx_func = func_string.split("(")[0].split(" ")[-1]
        if cxx_func not in self._known:
            logger.critical(f"Did not find {cxx_func!r} in doxygen tag file.")
            sys.exit(1)
        ix = 0
        if len(self._known[cxx_func]) > 1:
            # Disambiguate functions with same name
            args = re.findall(re.compile(r'(?<=\().*(?=\))'), func_string)
            if not args:
                known = '\n - '.join(
                    [""] + [ArgList.from_xml(xml_tag("arglist", xml)).short_str()
                            for xml in self._known[cxx_func]])
                # logger.debug(f"Known functions: are {known}")
                logger.critical(
                    f"Need argument list to disambiguate {func_string!r}. "
                    f"possible matches are:{known}")
                sys.exit(1)
            args = f"({args[0]})"
            ix = -1
            for i, xml in enumerate(self._known[cxx_func]):
                arglist = ArgList.from_xml(xml_tag("arglist", xml))
                if args[:-1] in arglist.short_str():
                    ix = i
                    break
            if ix < 0:
                logger.critical(
                    f"Unable to match {func_string!r} to known functions.")
                sys.exit(1)

        return TagInfo.from_xml(cxx_func, self._known[cxx_func][ix])


def tag_lookup(tag_info: TagInfo) -> TagDetails:
    """Retrieve tag details from doxygen tree."""
    xml_file = _xml_path / tag_info.anchorfile
    if not xml_file.exists():
        msg = (f"XML file does not exist at expected location: {xml_file}")
        logger.error(msg)
        return TagDetails()

    with xml_file.open() as fid:
        xml_details = fid.read()

    id_ = tag_info.id
    regex = re.compile(rf'<memberdef kind="function" id="{id_}"[\s\S]*?</memberdef>')
    matches = re.findall(regex, xml_details)

    if not matches:
        logger.error("No XML matches found for '%s'", tag_info.qualified_name)
        return TagDetails()
    if len(matches) != 1:
        logger.error("Inconclusive XML matches found for '%s'", tag_info.qualified_name)
        return TagDetails()

    def cleanup(entry: str) -> str:
        # Remove stray XML markup
        if entry.startswith("<para>"):
            entry = xml_tag("para", entry)
        if "<ref" in entry:
            regex = re.compile(r'<ref [\s\S]*?>')
            for ref in re.findall(regex, entry):
                entry = entry.replace(ref, "<ref>")
            entry = entry.replace("<ref>", "").replace("</ref>", "")
        return entry

    def resolve_parameteritem(par_map: str) -> List[Param]:
        # Resolve/flatten parameter list
        name_lines = xml_tag("parameternamelist", par_map).split("\n")
        regex = re.compile(r'(?<=<parametername)(.*?)(?=>)')
        names = []
        directions = []
        for name_line in name_lines:
            direction = re.findall(regex, name_line)[0]
            if "=" in direction:
                name_line = name_line.replace(direction, "")
                direction = direction.split("=")[1].strip('"')
            directions.append(direction)
            names.append(xml_tag("parametername", name_line))
        description = cleanup(xml_tag("parameterdescription", par_map))
        return [Param("", n, description, d) for n, d in zip(names, directions)]

    xml = matches[0]
    par_list = []
    par_block = xml_tag("parameterlist", xml, suffix='kind="param"')
    if par_block:
        for par_map in xml_tags("parameteritem", par_block):
            par_list.extend(resolve_parameteritem(par_map))

    def xml_attribute(attr: str, text: str, *, entry: str) -> str:
        """Extract XML attribute."""
        regex = re.compile(rf'(?<=<{attr} )(.*?)(?=/>)', flags=re.DOTALL)
        match = re.findall(regex, text)
        if not match:
            return ""  # not enough matches found
        entries = dict([tuple(_.split("=")) for _ in match[0].strip().split(" ")])
        return entries.get(entry, "").strip('"')

    return TagDetails(*tag_info,
                      cleanup(xml_tag("briefdescription", xml)),
                      xml_attribute("location", xml, entry="file"),
                      par_list)


def xml_tag(tag: str, text: str, suffix: str="", index=0) -> str:
    """Extract content enclosed between XML tags, optionally skipping matches."""
    if suffix:
        suffix = f" {suffix.strip()}"
    regex = re.compile(rf'(?<=<{tag}{suffix}>)(.*?)(?=</{tag}>)', flags=re.DOTALL)
    match = re.findall(regex, text)
    if index >= len(match):
        return ""  # not enough matches found
    return match[index].strip()


def xml_tags(tag: str, text: str, suffix: str="") -> List[str]:
    """Extract list of content enclosed by XML tags."""
    if suffix:
        suffix = f" {suffix.strip()}"
    regex = re.compile(rf'(?<=<{tag}{suffix}>)(.*?)(?=</{tag}>)',
                       flags=re.DOTALL|re.MULTILINE)
    matched = re.findall(regex, text)
    if not matched:
        blanks = text.split("\n")[-1].split("<")[0]
        msg = f"Could not extract {tag!r} from:\n{blanks}{text}\n"
        msg += f"using regex: {regex}"
        logger.error(msg)
        return []
    return matched
