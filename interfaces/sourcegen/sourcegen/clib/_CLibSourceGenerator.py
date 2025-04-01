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


_LOGGER = logging.getLogger()

class CLibSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating CLib."""

    _clib_bases: list[str] = None  #: list of bases provided via YAML configurations

    def __init__(self, out_dir: str, config: dict, templates: dict) -> None:
        self._out_dir = out_dir or None
        if self._out_dir is not None:
            self._out_dir = Path(out_dir)
            self._out_dir.mkdir(parents=True, exist_ok=True)
        self._config = Config.from_parsed(**config)
        self._templates = templates
        self._doxygen_tags = None

    @staticmethod
    def _javadoc_comment(block: str) -> str:
        """Build deblanked JavaDoc-style (C-style) comment block."""
        block = ["/**"] + block.strip().split("\n")
        block = "\n *  ".join(block).strip() + "\n */"
        return "\n".join([line.rstrip() for line in block.split("\n")])

    def _scaffold_annotation(self, c_func: CFunc, what: str) -> str:
        """Build annotation block via Jinja."""
        loader = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True)
        par_template = loader.from_string(self._templates["clib-param"])
        template = loader.from_string(self._templates["clib-comment"])

        def param(item: Param) -> str:
            ret = par_template.render(par=item)
            return f"{ret:<19} {item.description}"

        implements = what
        if isinstance(c_func.implements, CFunc):
            implements += f": {c_func.implements.short_declaration()}"
        block = template.render(
            brief=c_func.brief,
            params=[param(par) for par in c_func.arglist],
            returns=c_func.returns, implements=implements,
            relates=[f"{uu.base}::{uu.name}()" for uu in c_func.uses])
        return self._javadoc_comment(block)

    def _handle_crosswalk(
            self, what: str, crosswalk: dict, derived: dict[str, str]) -> str:
        """Crosswalk for object handle."""
        cabinet = None
        classes = list(self._config.includes.keys()) + list(derived.keys())
        for base in classes:
            ret_type = what.replace(f"<{base}>", "<T>")
            if ret_type in crosswalk:
                ret_type = crosswalk[ret_type]
                cabinet = base
                break
        if cabinet:
            # successful crosswalk with cabinet object
            return cabinet

        msg = f"Failed crosswalk for handle type {what!r} using {classes}."
        _LOGGER.critical(msg)
        sys.exit(1)

    def _ret_crosswalk(
            self, what: str, derived: list[str]) -> tuple[Param, list[Param]]:
        """Crosswalk for return type."""
        what = what.replace("virtual ", "")
        if what in self._config.ret_type_crosswalk:
            ret_type = self._config.ret_type_crosswalk[what]
            if ret_type == "void":
                returns = Param(
                    "int", "", "Zero for success or -1 for exception handling.")
                return returns, []
            if ret_type == "char*":
                # string expressions require special handling
                returns = Param(
                    "int", "", "Actual length of string including string-terminating "
                    "null byte, \\0, or -1 for exception handling.")
                buffer = [
                    Param("int", "bufLen", "Length of reserved array.", "in"),
                    Param(ret_type, "buf", "Returned string value.", "out")]
                return returns, buffer
            if ret_type.endswith("*"):
                # return type involves pointer to reserved buffer
                returns = Param(
                    "int", "",
                    "Actual length of value array or -1 for exception handling.")
                buffer = [
                    Param("int", "bufLen", "Length of reserved array.", "in"),
                    Param(ret_type, "buf", "Returned array value.", "out")]
                return returns, buffer
            if not any([what.startswith("shared_ptr"), ret_type.endswith("[]")]):
                # direct correspondence
                return Param(ret_type), []

            raise NotImplementedError(f"Crosswalk not implemented for {what!r}.")

        if "shared_ptr" in what:
            # check for crosswalk with object from includes
            handle = self._handle_crosswalk(
                what, self._config.ret_type_crosswalk, derived)
            returns = Param(
                "int", "",
                f"Handle to stored {handle} object or -1 for exception handling.")
            return returns, []

        msg = f"Failed crosswalk for return type {what!r}."
        _LOGGER.critical(msg)
        sys.exit(1)

    def _prop_crosswalk(self, par_list: list[Param]) -> list[Param]:
        """Crosswalk for argument type."""
        if not par_list:
            return []
        params = []
        for par in par_list:
            what = par.p_type
            if what in self._config.prop_type_crosswalk:
                if "vector<" in what:
                    params.append(
                        Param("int", f"{par.name}Len",
                              f"Length of vector reserved for {par.name}.", "in"))
                elif what.endswith("* const") or what.endswith("double*"):
                    direction = "in" if what.startswith("const") else "out"
                    params.append(
                        Param("int", f"{par.name}Len",
                              f"Length of array reserved for {par.name}.", direction))
                ret_type = self._config.prop_type_crosswalk[what]
                params.append(Param(ret_type, par.name, par.description, par.direction))
            elif "shared_ptr" in what:
                handle = self._handle_crosswalk(
                    what, self._config.prop_type_crosswalk, {})
                if "vector<" in what:
                    params.append(
                        Param("int", f"{par.name}Len",
                              f"Length of array reserved for {par.name}.", "in"))
                    description = f"Memory holding {handle} objects. "
                    description += par.description
                    params.append(Param("const int*", par.name, description.strip()))
                else:
                    description = f"Integer handle to {handle} object. "
                    description += par.description
                    params.append(
                        Param("int", par.name, description.strip(), par.direction))
            else:
                msg = f"Failed crosswalk for argument type {what!r}."
                _LOGGER.critical(msg)
                sys.exit(1)
        return params

    def _reverse_crosswalk(
            self, c_func: CFunc, base: str) -> tuple[dict[str, str], set[str]]:
        """Translate CLib arguments back to Jinja argument list."""
        handle = ""
        args = []
        lines = []
        buffer = []
        bases = set()

        def shared_object(cxx_type) -> str:
            """Extract object type from shared_ptr."""
            if "shared_ptr<" not in cxx_type:
                return None
            return cxx_type.split("<")[-1].split(">")[0]

        c_args = c_func.arglist
        cxx_member = c_func.implements
        if not cxx_member:
            if c_func.name.endswith("new"):
                # Default constructor
                cxx_func = CFunc("auto", f"make_shared<{base}>", ArgList([]))
            elif len(c_args) and "char*" in c_args[-1].p_type:
                cxx_func = CFunc("string", "dummy", ArgList([]), "", None, "", "base")
            else:
                cxx_func = CFunc("void", "dummy", ArgList([]), "", None, "", "base")
        elif isinstance(cxx_member, Param):
            if len(c_args) and "char*" in c_args[-1].p_type:
                cxx_func = CFunc("string", cxx_member.name, None,
                                 "", None, "", cxx_member.base)
            else:
                cxx_func = CFunc(cxx_member.p_type, cxx_member.name, None,
                                 "", None, "", cxx_member.base)
        else:
            cxx_func = cxx_member
        cxx_ix = 0
        check_array = False
        for c_ix, c_par in enumerate(c_func.arglist):
            c_name = c_par.name
            if isinstance(cxx_member, Param) or cxx_ix >= len(cxx_func.arglist):
                if c_ix == 0 and cxx_func.base and "len" not in c_name.lower():
                    handle = c_name
                    c_ix += 1
                if isinstance(cxx_member, Param) and cxx_member.direction == "out":
                    pass
                elif isinstance(cxx_member, Param):
                    break
                elif c_ix == len(c_args):
                    break

                # Handle output buffer and/or variable assignments
                cxx_type = cxx_func.ret_type
                if "string" in cxx_type:
                    buffer = ["string out",
                              f"copyString(out, {c_args[c_ix+1].name}, "
                              f"{c_args[c_ix].name});",
                              "int(out.size()) + 1"]  # include \0
                elif "vector" in cxx_type:
                    buffer = [f"{cxx_type} out",
                              "std::copy(out.begin(), out.end(), "
                              f"{c_args[c_ix+1].name});",
                              "int(out.size())"]
                elif "bool" in cxx_type:
                    buffer = [f"{cxx_type} out", "", "int(out)"]
                elif cxx_type in self._config.ret_type_crosswalk:
                    # can pass values directly
                    buffer = []
                else:
                    msg = (f"Scaffolding failed for {c_func.name!r}: reverse crosswalk "
                           f"not implemented for {cxx_type!r}:\n{c_func.declaration()}")
                    _LOGGER.critical(msg)
                    exit(1)
                break

            if isinstance(cxx_member, Param):
                cxx_arg = cxx_member
            else:
                cxx_arg = cxx_func.arglist[cxx_ix]
            if c_name != cxx_arg.name:
                # Encountered object handle or length indicator
                if c_name.endswith("Len"):
                    check_array = True
                elif c_ix == 0:
                    handle = c_name
                else:
                    msg = (f"Scaffolding failed for {c_func.name!r}: "
                           f"unexpected behavior for {c_name!r}.")
                    _LOGGER.critical(msg)
                    exit(1)
                continue

            cxx_type = cxx_arg.p_type
            if check_array:
                # Need to handle cross-walked parameter with length information
                c_prev = c_args[c_ix-1].name
                if "vector<shared_ptr" in cxx_type:
                    # Example: vector<shared_ptr<Domain1D>>
                    cxx_type = cxx_type.lstrip("const ").rstrip("&")
                    cxx_base = cxx_type.rstrip(">").split("<")[-1]
                    bases |= {cxx_base}
                    lines.extend([
                        f"{cxx_type} {c_name}_;",
                        f"for (int i = 0; i < {c_prev}; i++) {{",
                        f"    {c_name}_.push_back({cxx_base}Cabinet::at({c_name}[i]));",
                        "}",
                    ])
                    args.append(f"{c_name}_")
                elif "vector" in cxx_type:
                    # Example: vector<double> par_(par, par + parLen);
                    cxx_type = cxx_type.rstrip("&")
                    lines.append(
                        f"{cxx_type} {c_name}_({c_name}, {c_name} + {c_prev});")
                    args.append(f"{c_name}_")
                elif "*" in cxx_type:
                    # Can be passed directly; example: double *const
                    args.append(c_name)
                else:
                    msg = (f"Scaffolding failed for {c_func.name!r}: reverse "
                           f"crosswalk not implemented for {cxx_type!r}.")
                    _LOGGER.critical(msg)
                    exit(1)
                check_array = False
            elif "shared_ptr" in cxx_type:
                # Retrieve object from cabinet
                cxx_base = shared_object(cxx_type)
                args.append(f"{cxx_base}Cabinet::at({c_name})")
                if cxx_base != base:
                    bases |= {cxx_base}
            elif cxx_type == "bool":
                lines.append(f"bool {c_name}_ = ({c_name} != 0);")
                args.append(f"{c_name}_")
            else:
                # Regular parameter
                args.append(c_name)
            cxx_ix += 1

        # Obtain class and getter for managed objects
        shared = []
        checks = []
        for uu in c_func.uses:
            obj = shared_object(uu.ret_type)
            if obj:
                shared.append((obj, uu.name))
                bases |= {obj}
            else:
                checks.append(uu.name)

        # Ensure that all error codes are set correctly
        error = [-1, "ERR"]
        cxx_rtype = cxx_func.ret_type
        cxx_rbase = None
        if cxx_rtype.endswith(("int", "size_t")):
            error = ["ERR", "ERR"]
        elif cxx_rtype.endswith("double"):
            error = ["DERR", "DERR"]
        elif "shared_ptr" in cxx_rtype:
            cxx_rbase = shared_object(cxx_rtype)
            if cxx_rbase == base:
                buffer = ["auto obj", "", f"{cxx_rbase}Cabinet::index(*obj)"]
            else:
                buffer = ["auto obj", "", f"{cxx_rbase}Cabinet::index(*obj, {handle})"]
                bases |= {cxx_rbase}
            error = ["-2", "ERR"]
        elif cxx_rtype.endswith("void"):
            buffer = ["", "", "0"]

        if isinstance(cxx_member, Param) and cxx_member.direction == "in":
            c_name = c_func.arglist[-1].name
            if cxx_rtype.endswith("bool"):
                lines = [f"bool {c_name}_ = ({c_name} != 0);"]
                args.append(f"{c_name}_")
            else:
                args.append(c_name)

        ret = {
            "base": base, "handle": handle, "lines": lines, "buffer": buffer,
            "shared": shared, "checks": checks, "error": error, "cxx_rbase": cxx_rbase,
            "cxx_base": cxx_func.base, "cxx_name": cxx_func.name, "cxx_args": args,
            "cxx_implements": cxx_func.short_declaration(),
            "c_func": c_func.name, "c_args": [arg.name for arg in c_func.arglist],
        }
        return ret, bases

    def _scaffold_body(self, c_func: CFunc, recipe: Recipe) -> tuple[str, set[str]]:
        """Scaffold body of generic CLib function via Jinja."""
        loader = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True,
                             line_comment_prefix="##")
        args, bases = self._reverse_crosswalk(c_func, recipe.base)
        args["what"] = recipe.what

        if recipe.code:
            # override auto-generated code
            template = loader.from_string(self._templates["clib-custom-code"])
            args["lines"] = recipe.code.strip(" \n").split("\n")

        elif recipe.what == "noop":
            template = loader.from_string(self._templates["clib-noop"])

        elif recipe.what == "function":
            template = loader.from_string(self._templates["clib-function"])

        elif recipe.what == "variable-getter":
            template = loader.from_string(self._templates["clib-variable-getter"])

        elif recipe.what == "variable-setter":
            template = loader.from_string(self._templates["clib-variable-setter"])

        elif recipe.what == "constructor":
            template = loader.from_string(self._templates["clib-constructor"])

        elif recipe.what == "accessor":
            template = loader.from_string(self._templates["clib-accessor"])

        elif recipe.what == "destructor":
            template = loader.from_string(self._templates["clib-destructor"])

        elif recipe.what == "method":
            template = loader.from_string(self._templates["clib-method"])

        elif recipe.what == "getter":
            ret_type = c_func.implements.ret_type
            if "void" in ret_type or "vector" in ret_type:
                template = loader.from_string(self._templates["clib-array-getter"])
            else:
                template = loader.from_string(self._templates["clib-method"])

        elif recipe.what == "setter":
            if "*" in c_func.implements.arglist[0].p_type:
                template = loader.from_string(self._templates["clib-array-setter"])
            else:
                template = loader.from_string(self._templates["clib-method"])

        elif recipe.what == "reserved":
            args["cabinets"] = [kk for kk in self._clib_bases if kk]
            template = loader.from_string(
                self._templates[f"clib-reserved-{recipe.name}-cpp"])

        else:
            msg = f"{recipe.what!r} not implemented: {c_func.name!r}."
            _LOGGER.critical(msg)
            exit(1)

        body = template.render(**args)
        # remove blank lines left by line comments
        # see https://github.com/pallets/jinja/issues/204
        body = "\n".join(line for line in body.split("\n") if line.strip())
        return body, bases

    def _resolve_recipe(self, recipe: Recipe) -> CFunc:
        """Build CLib header from recipe and doxygen annotations."""
        def merge_params(
                implements: str, cxx_member: CFunc | Param
            ) -> tuple[list[Param], CFunc]:
            """Create preliminary CLib argument list."""
            obj_handle = []
            if "::" in implements:
                # If class method, add handle as first parameter
                what = implements.split("::")[0]
                obj_handle.append(
                    Param("int", "handle", f"Handle to queried {what} object."))
            if isinstance(cxx_member, Param):
                if recipe.what.endswith("setter"):
                    return obj_handle + [cxx_member], cxx_member
                return obj_handle, cxx_member

            if "(" not in implements:
                return obj_handle + cxx_member.arglist.params, cxx_member

            # Signature may skip C++ default parameters
            args_short = CFunc.from_str(implements).arglist
            if len(args_short) < len(cxx_member.arglist):
                cxx_arglist = ArgList(cxx_member.arglist[:len(args_short)])
                cxx_member = CFunc(cxx_member.ret_type, cxx_member.name,
                                   cxx_arglist, cxx_member.brief, cxx_member.implements,
                                   cxx_member.returns, cxx_member.base, cxx_member.uses)

            return obj_handle + cxx_member.arglist.params, cxx_member

        func_name = f"{recipe.prefix}_{recipe.name}"
        reserved = ["cabinetSize", "parentHandle",
                    "getCanteraError", "setLogWriter", "setLogCallback",
                    "clearStorage", "resetStorage"]
        if recipe.name in reserved:
            recipe.what = "reserved"
            loader = Environment(loader=BaseLoader)
            msg = f"   generating {func_name!r} -> {recipe.what}"
            _LOGGER.debug(msg)
            header = loader.from_string(
                self._templates[f"clib-reserved-{recipe.name}-h"]
                ).render(base=recipe.base, prefix=recipe.prefix)
            return CFunc.from_str(header, brief=recipe.brief)

        # Ensure that all functions/methods referenced in recipe are detected correctly
        bases = recipe.bases
        if not recipe.implements:
            recipe.implements = self._doxygen_tags.detect(recipe.name, bases)
        elif recipe.base and "::" not in recipe.implements:
            parts = list(recipe.implements.partition("("))
            parts[0] = self._doxygen_tags.detect(parts[0], bases)
            recipe.implements = "".join(parts)
        recipe.uses = [self._doxygen_tags.detect(uu.split("(")[0], bases, False)
                       for uu in recipe.uses]

        cxx_member = None
        ret_param = Param("void")
        args = []
        brief = ""

        if recipe.implements:
            cxx_member = self._doxygen_tags.cxx_member(
                recipe.implements, recipe.what.endswith("setter"))

            if cxx_member.base in recipe.derived:
                # Use alternative prefix for class specialization
                recipe.prefix = recipe.derived.get(cxx_member.base, recipe.prefix)
                func_name = f"{recipe.prefix}_{recipe.name}"

            msg = f"   generating {func_name!r} -> {recipe.implements}"
            _LOGGER.debug(msg)

            if isinstance(cxx_member, CFunc):
                # Convert C++ return type to format suitable for crosswalk:
                # Incompatible return parameters are buffered and appended to back
                ret_param, buffer_params = self._ret_crosswalk(
                    cxx_member.ret_type, recipe.derived)
                par_list, cxx_member = merge_params(recipe.implements, cxx_member)
                prop_params = self._prop_crosswalk(par_list)
                args = prop_params + buffer_params
                brief = cxx_member.brief
            elif recipe.what == "variable-setter":
                ret_param = Param("int")
                par_list, cxx_member = merge_params(recipe.implements, cxx_member)
                args = self._prop_crosswalk(par_list)
                brief = cxx_member.description
            else:
                # Variable getter
                prop_params, cxx_member = merge_params(recipe.implements, cxx_member)
                ret_param, buffer_params = self._ret_crosswalk(
                    cxx_member.p_type, recipe.derived)
                args = prop_params + buffer_params
                brief = cxx_member.description

        if recipe.what and cxx_member:
            # Recipe type and corresponding C++ function are known
            pass

        elif cxx_member and isinstance(cxx_member, Param):
            # Recipe represents a variable getter/setter
            recipe.what = "variable-getter"

        elif cxx_member:
            # Autodetection of CLib function purpose ("what")
            cxx_arglen = len(cxx_member.arglist)
            if not cxx_member.base:
                if (cxx_member.name.startswith("new") and
                    any(base in cxx_member.ret_type
                        for base in [recipe.base] + list(recipe.derived.keys()))):
                    recipe.what = "constructor"
                else:
                    recipe.what = "function"
            elif "void" not in cxx_member.ret_type and cxx_arglen == 0:
                recipe.what = "getter"
            elif "void" in cxx_member.ret_type and cxx_arglen == 1:
                p_type = cxx_member.arglist[0].p_type
                if cxx_member.name.startswith("get"):
                    recipe.what = "getter"
                elif "*" in p_type and not p_type.startswith("const"):
                    recipe.what = "getter"  # getter assigns to existing array
                else:
                    recipe.what = "setter"
            elif any(recipe.implements.startswith(base) for base in recipe.bases):
                recipe.what = "method"
            else:
                msg = f"Unable to auto-detect function type for recipe {recipe.name!r}."
                _LOGGER.critical(msg)
                exit(1)

        elif recipe.what == "noop":
            # No operation
            msg = f"   generating {func_name!r} -> no-operation"
            _LOGGER.debug(msg)
            brief = "No operation."
            ret_param = Param("int", "", "Always zero.")

        elif recipe.name == "new":
            # Default constructor
            recipe.what = "constructor"
            msg = f"   generating {func_name!r} -> default constructor"
            _LOGGER.debug(msg)
            brief= f"Instantiate {recipe.base} object using default constructor."
            ret_param = Param(
                "int", "", "Object handle if successful and -1 for exception handling.")

        elif recipe.name == "del":
            # Default destructor
            recipe.what = "destructor"
            msg = f"   generating {func_name!r} -> default destructor"
            _LOGGER.debug(msg)
            args = [Param("int", "handle", f"Handle to {recipe.base} object.")]
            brief= f"Delete {recipe.base} object."
            ret_param = Param(
                "int", "", "Zero for success and -1 for exception handling.")

        else:
            msg = f"Unable to resolve recipe type for {recipe.name!r}"
            _LOGGER.critical(msg)
            exit(1)

        if recipe.brief:
            brief = recipe.brief
        uses = [self._doxygen_tags.cxx_member(uu) for uu in recipe.uses]
        return CFunc(ret_param.p_type, func_name, ArgList(args), brief, cxx_member,
                     ret_param.description, None, uses)

    def _write_header(self, headers: HeaderFile) -> None:
        """Parse header specification and generate header file."""
        loader = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True)

        filename = headers.output_name(suffix="3.h")
        msg = f"  scaffolding {filename.name!r}"
        _LOGGER.info(msg)

        template = loader.from_string(self._templates["clib-definition"])
        declarations = []
        for c_func, recipe in zip(headers.funcs, headers.recipes):
            msg = f"   scaffolding {c_func.name!r} header"
            _LOGGER.debug(msg)
            declarations.append(
                template.render(
                    declaration=c_func.declaration(),
                    annotations=self._scaffold_annotation(c_func, recipe.what)))
        declarations = "\n\n".join(declarations)

        preamble = self._config.preambles.get(headers.base)

        guard = f"__{filename.name.upper().replace('.', '_')}__"
        t_file = Path(__file__).parent / "header_template.h.in"
        template = loader.from_string(t_file.read_text(encoding="utf-8"))
        output = template.render(
            name=filename.stem, guard=guard, preamble=preamble, prefix=headers.prefix,
            declarations=declarations, base=headers.base, docstring=headers.docstring)

        out = (Path(self._out_dir) /
               "include" / "cantera" / "clib_experimental" / filename.name)
        msg = f"  writing {filename.name!r}"
        _LOGGER.info(msg)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(output + "\n")

    def _write_implementation(self, headers: HeaderFile) -> None:
        """Parse header specification and generate implementation file."""
        loader = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True)

        filename = headers.output_name(suffix="3.cpp")
        msg = f"  scaffolding {filename.name!r}"
        _LOGGER.info(msg)

        template = loader.from_string(self._templates["clib-implementation"])
        implementations = []
        other = set()
        for c_func, recipe in zip(headers.funcs, headers.recipes):
            msg = f"   scaffolding {c_func.name!r} implementation"
            _LOGGER.debug(msg)
            body, bases = self._scaffold_body(c_func, recipe)
            implementations.append(
                template.render(declaration=c_func.declaration(),body=body))
            other |= bases
        implementations = "\n\n".join(implementations)
        str_utils = "copyString" in implementations

        if not headers.base:
            # main CLib file receives references to all cabinets
            other = [kk for kk in self._clib_bases if kk]
        includes = []
        for obj in [headers.base] + list(other):
            includes += self._config.includes[obj]

        t_file = Path(__file__).parent / "source_template.cpp.in"
        template = loader.from_string(t_file.read_text(encoding="utf-8"))
        output = template.render(
            name=filename.stem, implementations=implementations,
            prefix=headers.prefix, base=headers.base, docstring=headers.docstring,
            includes=includes, other=other, str_utils=str_utils)

        out = Path(self._out_dir) / "src" / "clib_experimental" / filename.name
        msg = f"  writing {filename.name!r}"
        _LOGGER.info(msg)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(output + "\n")

    def resolve_tags(self, headers_files: list[HeaderFile]) -> None:
        """Resolve recipe information based on doxygen tags."""
        def get_bases() -> tuple[list[str], list[str]]:
            bases = set()
            classes = set()
            for headers in headers_files:
                bases |= {headers.base}
                for recipe in headers.recipes:
                    classes |= set(recipe.bases)
            return sorted(bases), sorted(classes)

        self._clib_bases, classes = get_bases()
        self._doxygen_tags = TagFileParser(classes)

        for headers in headers_files:
            msg = f"  resolving recipes in {headers.path.name!r}:"
            _LOGGER.info(msg)
            c_funcs = []
            for recipe in headers.recipes:
                c_funcs.append(self._resolve_recipe(recipe))
            headers.funcs = c_funcs

    def generate_source(self, headers_files: list[HeaderFile]) -> None:
        """Generate output."""
        self.resolve_tags(headers_files)

        for headers in headers_files:
            self._write_header(headers)
            self._write_implementation(headers)
