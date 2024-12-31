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

    def __init__(self, out_dir: str, config: dict, templates: dict):
        self._out_dir = out_dir or None
        if self._out_dir is not None:
            self._out_dir = Path(out_dir)
            self._out_dir.mkdir(parents=True, exist_ok=True)
        self._config = Config.from_parsed(**config)  # typed config
        self._templates = templates
        self._doxygen_tags = None

    @staticmethod
    def _javadoc_comment(block):
        """Build deblanked JavaDoc-style (C-style) comment block."""
        block = ["/**"] + block.strip().split("\n")
        block = "\n * ".join(block).strip() + "\n */"
        return "\n".join([line.rstrip() for line in block.split('\n')])

    def _scaffold_annotation(self, c_func: CFunc, what: str) -> str:
        """Build annotation block via jinja."""
        loader = Environment(loader=BaseLoader)
        par_template = loader.from_string(self._templates["clib-param"])
        template = loader.from_string(self._templates["clib-comment"])

        def param(item: Param):
            ret = par_template.render(par=item)
            return f"{ret:<20} {item.description}"

        implements = what
        if isinstance(c_func.implements, CFunc):
            implements += f": {c_func.implements.short_declaration()}"
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
                elif what.endswith("*const"):
                    direction = "in" if what.startswith("const") else "out"
                    params.append(
                        Param("int", f"{par.name}Len",
                              f"Length of array reserved for {par.name}.", direction))
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

    @staticmethod
    def _reverse_crosswalk(
            c_args: list[Param], cxx_args: list[Param]) -> list[tuple[str, str]]:
        """Translate CLib arguments back to C++."""
        ret = []
        c_ix = 0
        for cxx_par in cxx_args:
            c_par = c_args[c_ix]
            if "shared_ptr" in cxx_par.p_type:
                # retrieve object from cabinet
                base = cxx_par.p_type.split("<")[-1].split(">")[0]
                ret.append((f"{base}Cabinet::at({c_par.name})", None))
            elif c_par.name.endswith("Len"):
                # need to create buffer variable
                ret.append((None, None))
                c_ix += 1
            else:
                # regular parameter
                ret.append((c_par.name, None))
            c_ix += 1
        return ret

    def _scaffold_getter(self, c_func: CFunc, recipe: Recipe) -> tuple[str, list[str]]:
        """Build getter via jinja."""
        loader = Environment(loader=BaseLoader)
        cxx_func = c_func.implements
        args = {
            "c_arg": [arg.name for arg in c_func.arglist],
            "cxx_base": recipe.base, "cxx_name": cxx_func.name,
            "cxx_implements": cxx_func.short_declaration(),
            "c_func": c_func.name, "size_fcn": recipe.size_fcn, "obj_base": None,
        }

        if "string" in cxx_func.ret_type:
            template = loader.from_string(self._templates["clib-string-getter"])
        elif "void" in cxx_func.ret_type:
            template = loader.from_string(self._templates["clib-array-getter"])
        elif "shared_ptr" in cxx_func.ret_type:
            args["obj_base"] = cxx_func.ret_type.split("<")[-1].split(">")[0]
            template = loader.from_string(self._templates["clib-object-getter"])
        elif "double" in c_func.ret_type:
            template = loader.from_string(self._templates["clib-double-getter"])
        elif "int" in c_func.ret_type:
            template = loader.from_string(self._templates["clib-int-getter"])
        else:
            logging.critical(f"Failed to scaffold getter: {args['cxx_implements']}")
            sys.exit(1)

        return template.render(**args), [args["obj_base"]]

    def _scaffold_setter(self, c_func: CFunc, recipe: Recipe) -> tuple[str, list[str]]:
        loader = Environment(loader=BaseLoader)
        cxx_func = c_func.implements
        args = {
            "c_arg": [arg.name for arg in c_func.arglist],
            "cxx_base": recipe.base, "cxx_name": cxx_func.name,
            "cxx_implements": cxx_func.short_declaration(),
            "c_func": c_func.name, "size_fcn": recipe.size_fcn, "obj_base": None,
        }

        p_type = cxx_func.arglist[0].p_type
        simple = ["string&", "int", "size_t", "double"]
        if "*" in p_type:
            template = loader.from_string(self._templates["clib-array-setter"])
        elif any(typ in p_type.split() for typ in simple):
            template = loader.from_string(self._templates["clib-simple-setter"])
        elif "shared_ptr" in p_type:
            args["obj_base"] = p_type.split("<")[-1].split(">")[0]
            template = loader.from_string(self._templates["clib-object-setter"])
        else:
            logging.critical(f"Failed to scaffold getter: {args['cxx_implements']}")
            sys.exit(1)

        return template.render(**args), [args["obj_base"]]

    def _scaffold_body(self, c_func: CFunc, recipe: Recipe) -> str:
        """Build function body via jinja."""

        loader = Environment(loader=BaseLoader)

        if not recipe.implements:
            return f"// {recipe.what}"

        if recipe.what == "getter":
            return self._scaffold_getter(c_func, recipe)[0]

        if recipe.what == "setter":
            return self._scaffold_setter(c_func, recipe)[0]

        cxx_func = c_func.implements
        return f"// {recipe.what}: {cxx_func.short_declaration()}"

    def _resolve_recipe(self, recipe: Recipe, quiet: bool=True) -> CFunc:
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

        if not recipe.implements:
            # autodetection of method
            for base in [recipe.base] + recipe.parents + recipe.derived:
                name = f"{base}::{recipe.name}"
                if self._doxygen_tags.exists(name):
                    recipe.implements = name
                    break

        cxx_func = None
        ret_param = Param("void")
        args = []
        brief = ""

        if recipe.implements:
            if not quiet:
                _logger.info(f"    generating {recipe.name!r} -> {recipe.implements}")
            cxx_func = self._doxygen_tags.cxx_func(recipe.implements, recipe.relates)

            # convert C++ return type to format suitable for crosswalk
            ret_param, buffer_params = self._ret_crosswalk(
                cxx_func.ret_type, recipe.derived)
            par_list = merge_params(recipe.implements, cxx_func)
            prop_params = self._prop_crosswalk(par_list)
            brief = cxx_func.brief
            args = prop_params + buffer_params

        if cxx_func and not recipe.what:
            # autodetection of CLib function purpose ("what")
            cxx_arglen = len(cxx_func.arglist)
            if any(base in cxx_func.ret_type
                   for base in [recipe.base] + recipe.derived) and \
                    cxx_func.name.startswith("new"):
                recipe.what = "constructor"
            elif "void" not in cxx_func.ret_type and cxx_arglen == 0:
                recipe.what = "getter"
            elif "void" in cxx_func.ret_type and cxx_arglen == 1:
                p_type = cxx_func.arglist[0].p_type
                if "*" in p_type and not p_type.startswith("const"):
                    recipe.what = "getter"  # getter assigns to existing array
                else:
                    recipe.what = "setter"
            elif any(recipe.implements.startswith(base)
                     for base in [recipe.base] + recipe.parents + recipe.derived):
                recipe.what = "method"
            else:
                recipe.what = "function"

        if recipe.what == "destructor":
            if not quiet:
                _logger.info(f"    generating {recipe.name!r} -> destructor")
            args = [Param("int", "handle", f"Handle to {recipe.base} object.")]
            brief= f"Delete {recipe.base} object."
            ret_param = Param(
                "int", "", "Zero for success and -1 for exception handling.")
            buffer_params = []

        func_name = f"{recipe.prefix}_{recipe.name}"
        return CFunc(ret_param.p_type, func_name, ArgList(args), brief, cxx_func,
                     ret_param.description)

    def _scaffold_header(self, header: HeaderFile) -> None:
        """Parse header specification and generate header file."""
        loader = Environment(loader=BaseLoader)

        template = loader.from_string(self._templates["clib-definition"])
        declarations = []
        for c_func, recipe in zip(header.funcs, header.recipes):
            declarations.append(
                template.render(
                    declaration=c_func.declaration(),
                    annotations=self._scaffold_annotation(c_func, recipe.what)))

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

    def _scaffold_source(self, header: HeaderFile) -> None:
        """Parse header specification and generate implementation file."""
        loader = Environment(loader=BaseLoader)

        template = loader.from_string(self._templates["clib-implementation"])
        implementations = []
        for c_func, recipe in zip(header.funcs, header.recipes):
            implementations.append(
                template.render(declaration=c_func.declaration(),
                                body=self._scaffold_body(c_func, recipe)))

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

    def resolve_tags(self, headers_files: list[HeaderFile], quiet: bool=True):
        """Resolve recipe information based on doxygen tags."""
        def get_bases() -> list[str]:
            bases = set()
            for headers in headers_files:
                for recipe in headers.recipes:
                    bases |= set([recipe.base] + recipe.uses + recipe.parents)
            return list(bases)
        self._doxygen_tags = TagFileParser(get_bases())

        for header in headers_files:
            if not quiet:
                _logger.info(f"  resolving recipes in {header.path.name!r}:")
            c_funcs = []
            for recipe in header.recipes:
                c_funcs.append(self._resolve_recipe(recipe, quiet=quiet))
            header.funcs = c_funcs

    def generate_source(self, headers_files: list[HeaderFile]):
        """Generate output."""
        self.resolve_tags(headers_files, quiet=False)

        for header in headers_files:
            self._scaffold_header(header)
            self._scaffold_source(header)
