"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from pathlib import Path
import logging

from jinja2 import Environment, BaseLoader

from ._Config import Config

from .._dataclasses import HeaderFile, Recipe #ArgList, Param, Func
from .._SourceGenerator import SourceGenerator
from .._TagFileParser import TagFileParser #, TagDetails, tag_lookup


_logger = logging.getLogger()

class CLibSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating CLib."""

    @staticmethod
    def _get_bases(headers_files: list[HeaderFile]) -> list[str]:
        bases = set()
        for headers in headers_files:
            for recipe in headers.recipes:
                bases |= set([recipe.base] + recipe.uses + recipe.parents)
        return list(bases)

    def build_declaration(self, recipe: Recipe) -> tuple[str, str, list[str]]:
        """Build strings containing declaration and annotation."""
        # def merge_params(implements, details: TagDetails) -> list[Param]:
        #     # If class method, add handle as first parameter
        #     args_merged = []
        #     if "::" in implements:
        #         what = implements.split("::")[0]
        #         args_merged.append(
        #             Param("int", "handle", f"Handle to queried {what} object."))
        #     # Merge parameters from signature, doxygen info and doxygen details
        #     args_used = ArgList.from_xml(details.arglist).params  # from doxygen
        #     if "(" in implements:
        #         args_short = Func.from_str(implements).arglist[:]  # from recipe
        #         args_used = args_used[:len(args_short)]
        #     args_annotated = details.parameterlist  # from documentation
        #     for arg in args_used:
        #         for desc in args_annotated:
        #             if arg.name == desc.name:
        #                 args_merged.append(
        #                     Param(arg.p_type, arg.name,
        #                           desc.description, desc.direction, arg.default))
        #                 break
        #         else:
        #             args_merged.append(arg)
        #     return args_merged

        if recipe.implements:
            tag_info = self._doxygen_tags.tag_info(recipe.implements)
            _logger.info(f"    recipe for {recipe.name!r}: {tag_info.anchorfile}")
            # details = tag_lookup(tag_info)

        #     # convert XML return type to format suitable for crosswalk
        #     ret_type = Param.from_xml(details.type).p_type
        #     ret_param, buffer_params, cabinets = self._ret_crosswalk(ret_type)
        #     par_list = merge_params(recipe.implements, details)
        #     prop_params, prop_cabinets = self._prop_crosswalk(par_list)
        #     cabinets += prop_cabinets
        #     args = prop_params + buffer_params

        elif recipe.what == "destructor":
            _logger.info(f"    recipe for {recipe.name!r}: destructor")
        #     args = [Param("int", "handle", f"Handle to {recipe.base} object.")]
        #     details = TagDetails(
        #         "", "", "", "", "", "", f"Delete {recipe.base} object.", "", args)
        #     ret_param = Param(
        #         "int", "", "Zero for success and -1 for exception handling.")
        #     annotations = f"//! {details.briefdescription}"
        #     buffer_params = []
        #     cabinets = [recipe.base]

        else:
            _logger.critical("Unable to build declaration for '%s' with type '%s'.",
                            recipe.name, recipe.what)
            sys.exit(1)

        # c_func = Func(ret_param.p_type, recipe.name, ArgList(args))
        # declaration = c_func.declaration()
        # annotations = self._build_annotation(details, ret_param, args, recipe.relates)

        # return declaration, annotations, cabinets
        return "", "", []


    def _parse_header(self, header: HeaderFile):
        """Parse header file and generate output."""
        loader = Environment(loader=BaseLoader)

        template = loader.from_string(self._templates["clib-definition"])
        declarations = []
        for recipe in header.recipes:
            declaration, annotations, _ = self.build_declaration(recipe)
        #     declarations.append(
        #         template.render(declaration=declaration, annotations=annotations))

        filename = header.output_name(suffix=".h", auto="3")
        guard = f"__{filename.name.upper().replace('.', '_')}__"
        template = loader.from_string(self._templates["clib-header-file"])
        output = template.render(
            name=filename.stem, guard=guard, header_entries=declarations)

        if self._out_dir:
            out = Path(self._out_dir) / "include" / filename.name
            _logger.info(f"  writing {filename!r}")
            if not out.parent.exists():
                out.parent.mkdir(parents=True)
            with open(out, "wt", encoding="utf-8") as stream:
                stream.write(output)
                stream.write("\n")
        else:
            print(output)

    def __init__(self, out_dir: str, config: dict, templates: dict):
        self._out_dir = out_dir or None
        if self._out_dir is not None:
            self._out_dir = Path(out_dir)
            self._out_dir.mkdir(parents=True, exist_ok=True)
        self._config = Config.from_parsed(**config)  # typed config
        self._templates = templates

    def generate_source(self, headers_files: list[HeaderFile]):
        """Generate output."""
        self._doxygen_tags = TagFileParser(self._get_bases(headers_files))

        for header in headers_files:
            _logger.info(f"  parsing recipes in {header.path.name!r}:")
            self._parse_header(header)
