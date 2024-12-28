"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from pathlib import Path
import logging

from jinja2 import Environment, BaseLoader

from ._Config import Config

from .._dataclasses import HeaderFile, Recipe, Param, ArgList, Func
from .._SourceGenerator import SourceGenerator
from .._TagFileParser import TagFileParser, TagDetails, tag_lookup


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

    @staticmethod
    def _javadoc_comment(block):
        """Build deblanked JavaDoc-style (C-style) comment block."""
        block = ["/**"] + block.strip().split("\n")
        block = "\n * ".join(block).strip() + "\n */"
        return "\n".join([line.rstrip() for line in block.split('\n')])

    def _build_annotation(self, details: TagDetails, ret: Param,
                          params: list[Param], relates: list=[]) -> str:
        """Build annotation block."""
        def param(item: Param):
            ret = par_template.render(par=item)
            return f"{ret:<20} {item.description or 'Undocumented.'}"

        loader = Environment(loader=BaseLoader)
        par_template = loader.from_string(self._templates["clib-param"])
        template = loader.from_string(self._templates["clib-comment"])
        implements = ""
        if details.qualified_name:
            arglist = ArgList.from_xml(details.arglist)
            implements = f"{details.qualified_name}{arglist.short_str()}"

        block = template.render(
            briefdescription=details.briefdescription,
            params=[param(par) for par in params],
            returns=ret.description, implements=implements,
            relates=relates)
        return self._javadoc_comment(block)

    def _handle_crosswalk(self, what: str, crosswalk: dict, derived: list[str]) -> str:
        """Crosswalk for object handle."""
        cabinet = None
        classes = self._config.cabinets + derived
        for base in classes:
            ret_type = what.replace(f"<{base}>", "<T>")
            if ret_type in crosswalk:
                ret_type = crosswalk[ret_type]
                cabinet = base
                break
        if cabinet:
            # successful crosswalk with cabinet object
            return cabinet

        logging.critical(f"Failed crosswalk for handle type {what!r} using {classes}.")
        sys.exit(1)

    def _ret_crosswalk(self, what: str, derived: list[str]) -> tuple[Param, list[Param], list[str]]:
        """Crosswalk for return type."""
        what = what.replace("virtual ", "")
        if what in self._config.ret_type_crosswalk:
            ret_type = self._config.ret_type_crosswalk[what]
            if ret_type == "char*":
                # string expressions require special handling
                returns = Param(
                    "int", "",
                    "Actual length of string or -1 for exception handling.")
                buffer = [
                    Param("int", "lenBuf", "Length of reserved array.", "in"),
                    Param(ret_type, "charBuf", "Returned string value.", "out")]
                return returns, buffer, []
            if ret_type == "void":
                returns = Param(
                    "int", "",
                    "Zero for success or -1 for exception handling.")
                return returns, [], []
            if not any([what.startswith("shared_ptr"), ret_type.endswith("[]")]):
                # direct correspondence
                return Param(ret_type), [], []

            # all other types require reserved buffers
            what = "vector" if what.startswith("vector") else "array"
            returns = Param(
                "int", "", f"Actual length of {what} or -1 for exception handling.")
            buffer = [
                Param("int", "lenBuf", "Length of reserved array.", "in"),
                Param(ret_type, "valueBuf", f"Returned {what} value.", "out")]
            return returns, buffer, []

        if "shared_ptr" in what:
            # check for crosswalk with object from cabinets
            handle = self._handle_crosswalk(
                what, self._config.ret_type_crosswalk, derived)
            returns = Param(
                "int", "",
                f"Handle to stored {handle} object or -1 for exception handling.")
            return returns, [], [handle]

        logging.critical(f"Failed crosswalk for return type {what!r}.")
        sys.exit(1)

    def _prop_crosswalk(self, par_list: list[Param]) -> tuple[list[Param], list[str]]:
        """Crosswalk for argument type."""
        if not par_list:
            return [], []
        params = []
        cabinets = set()
        for par in par_list:
            what = par.p_type
            if what in self._config.prop_type_crosswalk:
                if "vector" in what:
                    params.append(
                        Param("int", f"{par.name}Len",
                              f"Length of vector reserved for {par.name}.", "in"))
                ret_type = self._config.prop_type_crosswalk[what]
                params.append(Param(ret_type, par.name, par.description, par.direction))
            elif "shared_ptr" in what:
                handle = self._handle_crosswalk(what, self._config.prop_type_crosswalk, [])
                cabinets |= {handle}
                params.append(Param("int", par.name, par.description, par.direction))
            else:
                logging.critical(f"Failed crosswalk for argument type {what!r}.")
                sys.exit(1)
        return params, list[cabinets]

    def build_declaration(self, recipe: Recipe) -> tuple[str, str, list[str]]:
        """Build strings containing declaration and annotation."""
        def merge_params(implements, details: TagDetails) -> list[Param]:
            # If class method, add handle as first parameter
            args_merged = []
            if "::" in implements:
                what = implements.split("::")[0]
                args_merged.append(
                    Param("int", "handle", f"Handle to queried {what} object."))
            # Merge parameters from signature, doxygen info and doxygen details
            args_used = ArgList.from_xml(details.arglist).params  # from doxygen
            if "(" in implements:
                args_short = Func.from_str(implements).arglist[:]  # from recipe
                args_used = args_used[:len(args_short)]
            args_annotated = details.parameterlist  # from documentation
            for arg in args_used:
                for desc in args_annotated:
                    if arg.name == desc.name:
                        args_merged.append(
                            Param(arg.p_type, arg.name,
                                  desc.description, desc.direction, arg.default))
                        break
                else:
                    args_merged.append(arg)
            return args_merged

        if recipe.implements:
            _logger.info(f"    generating {recipe.name!r} -> {recipe.implements}")
            tag_info = self._doxygen_tags.tag_info(recipe.implements)
            details = tag_lookup(tag_info)

            # convert XML return type to format suitable for crosswalk
            ret_type = Param.from_xml(details.type).p_type
            ret_param, buffer_params, cabinets = self._ret_crosswalk(
                ret_type, recipe.derived)
            par_list = merge_params(recipe.implements, details)
            prop_params, prop_cabinets = self._prop_crosswalk(par_list)
            cabinets += prop_cabinets
            args = prop_params + buffer_params

        elif recipe.what == "destructor":
            _logger.info(f"    generating {recipe.name!r} -> destructor")
            args = [Param("int", "handle", f"Handle to {recipe.base} object.")]
            details = TagDetails(
                "", "", "", "", "", "", "", f"Delete {recipe.base} object.", args)
            ret_param = Param(
                "int", "", "Zero for success and -1 for exception handling.")
            annotations = f"//! {details.briefdescription}"
            buffer_params = []
            cabinets = [recipe.base]

        else:
            _logger.critical(f"Unable to build declaration for {recipe.name!r} "
                             f"with type {recipe.what!r}.")
            sys.exit(1)

        c_func = Func(ret_param.p_type, recipe.name, ArgList(args))
        declaration = c_func.declaration()
        annotations = self._build_annotation(details, ret_param, args, recipe.relates)

        return declaration, annotations, cabinets

    def _parse_header(self, header: HeaderFile):
        """Parse header file and generate output."""
        loader = Environment(loader=BaseLoader)

        template = loader.from_string(self._templates["clib-definition"])
        declarations = []
        for recipe in header.recipes:
            declaration, annotations, _ = self.build_declaration(recipe)
            declarations.append(
                template.render(declaration=declaration, annotations=annotations))

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
