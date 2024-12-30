"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from pathlib import Path
import logging

from jinja2 import Environment, BaseLoader

from ._Config import Config

from .._dataclasses import HeaderFile, Param, ArgList, CFunc, Recipe
from .._SourceGenerator import SourceGenerator
from .._TagFileParser import TagFileParser


_logger = logging.getLogger()

class CLibSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating CLib."""

    @staticmethod
    def _javadoc_comment(block):
        """Build deblanked JavaDoc-style (C-style) comment block."""
        block = ["/**"] + block.strip().split("\n")
        block = "\n * ".join(block).strip() + "\n */"
        return "\n".join([line.rstrip() for line in block.split('\n')])

    def _build_annotation(self, c_func: CFunc) -> str:
        """Build annotation block via jinja."""
        loader = Environment(loader=BaseLoader)
        par_template = loader.from_string(self._templates["clib-param"])
        template = loader.from_string(self._templates["clib-comment"])

        def param(item: Param):
            ret = par_template.render(par=item)
            return f"{ret:<20} {item.description}"

        implements = ""
        if isinstance(c_func.implements, CFunc):
            implements = c_func.implements.short_declaration()
        block = template.render(
            briefdescription=c_func.brief,
            params=[param(par) for par in c_func.arglist],
            returns=c_func.returns, implements=implements, relates=c_func.relates)
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

    def _ret_crosswalk(self, what: str, derived: list[str]) -> tuple[Param, list[Param]]:
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
                return returns, buffer
            if ret_type == "void":
                returns = Param(
                    "int", "", "Zero for success or -1 for exception handling.")
                return returns, []
            if not any([what.startswith("shared_ptr"), ret_type.endswith("[]")]):
                # direct correspondence
                return Param(ret_type), []

            # all other types require reserved buffers
            what = "vector" if what.startswith("vector") else "array"
            returns = Param(
                "int", "", f"Actual length of {what} or -1 for exception handling.")
            buffer = [
                Param("int", "lenBuf", "Length of reserved array.", "in"),
                Param(ret_type, "valueBuf", f"Returned {what} value.", "out")]
            return returns, buffer

        if "shared_ptr" in what:
            # check for crosswalk with object from cabinets
            handle = self._handle_crosswalk(
                what, self._config.ret_type_crosswalk, derived)
            returns = Param(
                "int", "",
                f"Handle to stored {handle} object or -1 for exception handling.")
            return returns, []

        logging.critical(f"Failed crosswalk for return type {what!r}.")
        sys.exit(1)

    def _prop_crosswalk(self, par_list: list[Param]) -> list[Param]:
        """Crosswalk for argument type."""
        if not par_list:
            return [], []
        params = []
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
                handle = self._handle_crosswalk(
                    what, self._config.prop_type_crosswalk, [])
                description = f"Integer handle to {handle} object. {par.description}"
                params.append(Param("int", par.name, description, par.direction))
            else:
                logging.critical(f"Failed crosswalk for argument type {what!r}.")
                sys.exit(1)
        return params

    def clib_header(self, recipe: Recipe) -> tuple[str, str, list[str]]:
        """Build CLib header from recipe and doxygen annotations."""
        def merge_params(implements, cxx_func: CFunc) -> list[Param]:
            # If class method, add handle as first parameter
            args_merged = []
            if "::" in implements:
                what = implements.split("::")[0]
                args_merged.append(
                    Param("int", "handle", f"Handle to queried {what} object."))
            # Combine signature and annotated C++ function
            args_annotated = cxx_func.arglist.params
            if "(" in implements:
                args_short = CFunc.from_str(implements).arglist[:]  # from recipe
                args_annotated = args_annotated[:len(args_short)]
            return args_merged + args_annotated

        if recipe.implements:
            _logger.info(f"    generating {recipe.name!r} -> {recipe.implements}")
            cxx_func = self._doxygen_tags.cxx_func(recipe.implements, recipe.relates)

            # convert XML return type to format suitable for crosswalk
            ret_param, buffer_params = self._ret_crosswalk(
                cxx_func.ret_type, recipe.derived)
            par_list = merge_params(recipe.implements, cxx_func)
            prop_params = self._prop_crosswalk(par_list)
            brief = cxx_func.brief
            args = prop_params + buffer_params

        elif recipe.what == "destructor":
            _logger.info(f"    generating {recipe.name!r} -> destructor")
            args = [Param("int", "handle", f"Handle to {recipe.base} object.")]
            brief= f"Delete {recipe.base} object."
            cxx_func = None
            ret_param = Param(
                "int", "", "Zero for success and -1 for exception handling.")
            buffer_params = []

        else:
            _logger.critical(f"Unable to build declaration for {recipe.name!r} "
                             f"with type {recipe.what!r}.")
            sys.exit(1)

        return CFunc(ret_param.p_type, recipe.name, ArgList(args), brief, cxx_func,
                     ret_param.description)

    def build_header(self, header: HeaderFile) -> None:
        """Parse header specification and generate header file."""
        loader = Environment(loader=BaseLoader)

        template = loader.from_string(self._templates["clib-definition"])
        declarations = []
        for recipe in header.recipes:
            c_func = self.clib_header(recipe)
            declarations.append(
                template.render(declaration=c_func.declaration(),
                                annotations=self._build_annotation(c_func)))

        filename = header.output_name(suffix=".h", auto="3")
        guard = f"__{filename.name.upper().replace('.', '_')}__"
        template = loader.from_string(self._templates["clib-header-file"])
        output = template.render(
            name=filename.stem, guard=guard, header_entries=declarations)

        if self._out_dir:
            out = Path(self._out_dir) / "include" / filename.name
            _logger.info(f"  writing {filename.name!r}")
            if not out.parent.exists():
                out.parent.mkdir(parents=True, exist_ok=True)
            with open(out, "wt", encoding="utf-8") as stream:
                stream.write(output)
                stream.write("\n")
        else:
            print(output)

    def build_source(self, header: HeaderFile) -> None:
        """Parse header specification and generate implementation file."""
        loader = Environment(loader=BaseLoader)

        template = loader.from_string(self._templates["clib-implementation"])
        implementations = []
        for recipe in header.recipes:
            c_func = self.clib_header(recipe)
            implementations.append(
                template.render(declaration=c_func.declaration()))

        filename = header.output_name(suffix=".cpp", auto="3")
        template = loader.from_string(self._templates["clib-source-file"])
        output = template.render(
            name=filename.stem, source_entries=implementations)

        if self._out_dir:
            out = Path(self._out_dir) / "src" / filename.name
            _logger.info(f"  writing {filename.name!r}")
            if not out.parent.exists():
                out.parent.mkdir(parents=True, exist_ok=True)
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
        self._doxygen_tags = None

    def parse_tags(self, headers_files: list[HeaderFile]):
        """Parse doxygen tags."""
        def get_bases() -> list[str]:
            bases = set()
            for headers in headers_files:
                for recipe in headers.recipes:
                    bases |= set([recipe.base] + recipe.uses + recipe.parents)
            return list(bases)
        self._doxygen_tags = TagFileParser(get_bases())

    def generate_source(self, headers_files: list[HeaderFile]):
        """Generate output."""
        self.parse_tags(headers_files)

        for header in headers_files:
            _logger.info(f"  parsing recipes in {header.path.name!r}:")
            self.build_header(header)
            self.build_source(header)
