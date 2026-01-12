"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging
from dataclasses import dataclass

from jinja2 import Environment, BaseLoader

from ..dataclasses import HeaderFile, Param, ArgList, Func, Recipe
from ..generator import SourceGenerator
from .._helpers import with_unpack_iter, escape_token


_LOGGER = logging.getLogger()


@dataclass(frozen=True)
@with_unpack_iter
class Config:
    """Provides configuration info for the CLibSourceGenerator class"""

    ret_type_crosswalk: dict[str, str]  #: Return type crosswalks

    par_type_crosswalk: dict[str, str]  #: Parameter type crosswalks

    preambles: dict[str, str]  #: Preamble text for each header file

    includes: dict[str, list[str]]  #: Include directives for each implementation file


class CLibSourceGenerator(SourceGenerator):
    """The SourceGenerator for generating CLib."""

    _clib_bases: list[str] = None  #: list of bases provided via YAML configurations

    def __init__(self, out_dir: str, config: dict, templates: dict) -> None:
        self._out_dir = out_dir or None
        if self._out_dir is not None:
            self._out_dir = Path(out_dir)
            self._out_dir.mkdir(parents=True, exist_ok=True)
        self._config = Config(**config)
        self._templates = templates
        self._doxygen_tags = None

    @staticmethod
    def _javadoc_comment(block: str) -> str:
        """Build deblanked JavaDoc-style (C-style) comment block."""
        block = ["/**"] + block.strip().split("\n")
        block = "\n *  ".join(block).strip() + "\n */"
        return "\n".join([line.rstrip() for line in block.split("\n")])

    def _scaffold_annotation(self, c_func: Func, what: str) -> str:
        """Build annotation block via Jinja."""
        loader = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True)
        par_template = loader.from_string(self._templates["clib-param"])
        template = loader.from_string(self._templates["clib-comment"])

        def param(item: Param) -> str:
            ret = par_template.render(par=item)
            return f"{ret:<19} {item.description}"

        if isinstance(c_func.wraps, Func):
            wraps = c_func.wraps.short_declaration()
        elif isinstance(c_func.wraps, Param):
            wraps = c_func.wraps.long_str()
        elif isinstance(c_func.wraps, str):
            wraps = "custom code"
        else:
            wraps = "undefined"
        uses = [ uu.short_declaration() for uu in c_func.uses]
        block = template.render(
            brief=c_func.brief,
            params=[param(par) for par in c_func.arglist],
            returns=c_func.returns, what=what, wraps=wraps, uses=uses,
            deprecated=c_func.deprecated)
        return self._javadoc_comment(block)

    @staticmethod
    def _shared_object(cxx_type) -> str:
        """Extract object type from shared_ptr."""
        if "shared_ptr<" not in cxx_type:
            return None
        return cxx_type.split("<")[-1].split(">")[0]

    @staticmethod
    def _critical(c_func: Func, cxx_func: Func, msg) -> None:
        msg = (f"Scaffolding failed for {c_func.name!r}: {msg}:\n"
                f"C: {c_func.declaration()}\n")
        if isinstance(cxx_func.wraps, str):
            msg += f"C++: {cxx_func.wraps} (wrapped as {cxx_func.declaration()})"
        else:
            msg += f"C++: {cxx_func.declaration()}"
        _LOGGER.critical(msg)
        exit(1)

    def _ret_crosswalk(self, c_func: Func, cxx_func: Func
                       ) -> tuple[str, list, str | None, str | None]:
        """Crosswalk for C++ return type."""
        c_args = c_func.arglist
        c_type = c_func.ret_type
        handle = ""
        buffer = []
        after = None
        cxx_rbase = None

        cxx_type = cxx_func.ret_type
        if cxx_type.startswith("virtual "):
            cxx_type = cxx_type.replace("virtual ", "")
        ret_key = cxx_type.removeprefix("const ").removesuffix(" const").rstrip("&")

        if c_args and cxx_func.base:
            handle = c_args[0].name

        # check C++ type information
        shared_base = self._shared_object(ret_key)
        cxx_template = ret_key.replace(shared_base, "T") if shared_base else None
        if cxx_template:
            ret_key = cxx_template
            if ret_key not in self._config.par_type_crosswalk:
                self._critical(
                    c_func, cxx_func,
                    f"return crosswalk not listed for C++ template {ret_key!r}")
        elif cxx_type in ["auto"] or cxx_func.wraps == "custom code":
            pass
        elif ret_key not in self._config.ret_type_crosswalk:
            self._critical(c_func, cxx_func,
                           f"return crosswalk not listed for C++ type {ret_key!r}")

        # check C type information
        if cxx_type == "auto" or cxx_func.wraps == "custom code":
            pass
        elif c_args and c_args[-1].name == "buf":
            # return buffer
            if c_type != "int32_t":
                self._critical(c_func, cxx_func,
                               "method with buffered return should return 'int32_t'")
            c_type = c_args[-1].p_type
            c_key = c_type.removeprefix("const ")
            if c_key not in self._config.ret_type_crosswalk[ret_key]:  # 1
                self._critical(
                    c_func, cxx_func,
                    f"buffered return crosswalk not listed for C type {c_key!r}")
        else:
            c_key = c_type.removeprefix("const ")
            if c_key not in self._config.ret_type_crosswalk[ret_key]:  # 2
                self._critical(c_func, cxx_func,
                               f"return crosswalk not listed for C type {c_key!r}")

        if cxx_type == "bool":
            buffer = [f"{cxx_type} out", "int(out)"]
        elif cxx_type.endswith("void"):
            buffer = ["", "0"]
        elif "string" in cxx_type:
            buffer = [
                "string out",
                "int(out.size()) + 1",  # include \0
            ]
            after = f"copyString(out, {c_args[-1].name}, {c_args[-2].name});"
        elif "shared_ptr" in cxx_type:
            cxx_rbase = self._shared_object(cxx_type)
        elif "vector" in cxx_type:
            buffer = [
                f"{cxx_type} out",
                "int(out.size())",
            ]
            after = f"std::copy(out.begin(), out.end(), {c_args[-1].name});"
        elif cxx_type.endswith("*"):
            buffer = [
                "auto out",
                c_args[-2].name,
            ]
            after = f"std::copy(out, out + {c_args[-2].name}, {c_args[-1].name});"

        return handle, buffer, after, cxx_rbase

    def _par_crosswalk(self, c_func: Func, cxx_func: Func, base: str
                       ) -> tuple[list, list, list, set]:
        """Crosswalk for C++ return type."""
        args = []
        before = []
        after = []
        bases = set()

        c_args = c_func.arglist
        cxx_args = cxx_func.arglist

        c_ix = 0
        for cxx_arg in cxx_args:
            while c_args[c_ix].name != cxx_arg.name:
                c_ix += 1
            check_array = False
            if c_ix > 0 and c_args[c_ix-1].name.endswith("Len"):
                c_type = c_args[c_ix-1].p_type
                if c_type != "int32_t":
                    self._critical(c_func, cxx_func,
                                   f"invalid variable setter return type {c_type!r}")
                check_array = True

            c_name = c_args[c_ix].name
            c_type = c_args[c_ix].p_type
            cxx_type = cxx_arg.p_type.removesuffix(" const")

            par_key = cxx_type.removeprefix("const ").removesuffix(" const").rstrip("&")
            shared_base = self._shared_object(par_key)
            cxx_template = par_key.replace(shared_base, "T") if shared_base else None
            if cxx_template:
                par_key = cxx_template
                if par_key not in self._config.par_type_crosswalk:
                    self._critical(
                        c_func, cxx_func,
                        f"crosswalk not listed for C++ template {par_key!r}")
                par_key = cxx_template
            elif par_key not in self._config.par_type_crosswalk:
                self._critical(c_func, cxx_func,
                               f"crosswalk not listed for C++ type {par_key!r}")

            c_key = c_type.removeprefix("const ")
            if c_key not in self._config.par_type_crosswalk[par_key]:
                self._critical(c_func, cxx_func,
                               f"crosswalk not listed for C type {c_key!r}")

            if check_array:
                # Need to handle crosswalked parameter with length information
                c_prev = c_args[c_ix-1].name
                if "vector<shared_ptr" in cxx_type:
                    # Example: vector<shared_ptr<Domain1D>>
                    cxx_type = cxx_type.removeprefix("const ").rstrip("&")
                    cxx_base = cxx_type.rstrip(">").split("<")[-1]
                    bases.add(cxx_base)
                    before.extend([
                        f"{cxx_type} {c_name}_;",
                        f"for (int i = 0; i < {c_prev}; i++) {{",
                        f"    {c_name}_.push_back({cxx_base}Cabinet::at({c_name}[i]));",
                        "}",
                    ])
                    args.append(f"{c_name}_")
                elif "vector<string" in cxx_type:
                    # requires conversion to character arrays or YAML serialization
                    self._critical(c_func, cxx_func,
                                   f"crosswalk not implemented for C++ {cxx_type!r}.")
                elif "vector" in cxx_type:
                    # Example: vector<double> par_(par, par + parLen);
                    cxx_type = cxx_type.rstrip("&")
                    if cxx_arg.direction == "out":
                        before.append(f"{cxx_type} {c_name}_({c_prev});")
                    else:
                        before.append(
                            f"{cxx_type} {c_name}_({c_name}, {c_name} + {c_prev});")
                    if not cxx_type.startswith("const") or cxx_arg.direction == "out":
                        after.append(
                            f"std::copy({c_name}_.begin(), {c_name}_.end(), {c_name});")
                    args.append(f"{c_name}_")
                elif "span<" in cxx_type:
                    # Example: span<double> par_(par, par + parLen);
                    before.append(f"{cxx_type} {c_name}_({c_name}, {c_prev});")
                    args.append(f"{c_name}_")
                elif "*" in cxx_type:
                    # Can be passed directly; example: double *const
                    args.append(c_name)
                else:
                    self._critical(c_func, cxx_func,
                                   f"crosswalk not implemented for C++ {cxx_type!r}.")
            elif "shared_ptr" in cxx_type:
                # Retrieve object from cabinet
                args.append(f"{shared_base}Cabinet::at({c_name})")
                if shared_base != base:
                    bases.add(shared_base)
            elif cxx_type == "bool":
                before.append(f"bool {c_name}_ = ({c_name} != 0);")
                args.append(f"{c_name}_")
            else:
                # Regular parameter
                args.append(c_name)

        return args, before, after, bases

    def _crosswalk(self, c_func: Func, base: str) -> tuple[dict[str, str], set[str]]:
        """Translate CLib arguments back to Jinja argument list."""
        c_args = c_func.arglist
        cxx_wraps = None
        cxx_member = c_func.wraps
        if not cxx_member:
            if c_func.name.endswith("new"):
                cxx_wraps = f"default constructor for {base}"
                cxx_func = Func("auto", f"make_shared<{base}>", ArgList([]),
                                "", cxx_wraps, "", base)
            elif len(c_args) and "char*" in c_args[-1].p_type:
                cxx_wraps = "reserved function/method"
                cxx_func = Func("string", "dummy", ArgList([]),
                                "", cxx_wraps, "", base)
            else:
                cxx_wraps = f"reserved function/method: {c_func.name}"
                cxx_func = Func("void", "dummy", ArgList([]),
                                "", cxx_wraps, "", base)
        elif isinstance(cxx_member, Param):
            if cxx_member.direction == "in":
                cxx_wraps = f"variable setter: {cxx_member.name}"
                cxx_func = Func("int", cxx_member.name, ArgList([cxx_member]),
                                "", cxx_wraps, "", cxx_member.base)
            else:  # cxx_member.direction == "out"
                cxx_wraps = f"variable getter: {cxx_member.name}"
                if len(c_args) and "char*" in c_args[-1].p_type:
                    cxx_func = Func("string", cxx_member.name, ArgList([]),
                                    "", cxx_wraps, "", cxx_member.base)
                else:
                    cxx_func = Func(cxx_member.p_type, cxx_member.name, ArgList([]),
                                    "", cxx_wraps, "", cxx_member.base)
        elif isinstance(cxx_member, str):
            cxx_wraps = "custom code"
            cxx_func = Func("void", "dummy", ArgList([]), "", cxx_wraps, "", base)
        else:
            cxx_wraps = cxx_member.short_declaration()
            cxx_func = cxx_member

        handle, buffer, rafter, rbase = self._ret_crosswalk(c_func, cxx_func)
        args, before, after, bases = self._par_crosswalk(c_func, cxx_func, base)
        if rafter:
            after.append(rafter)
        if rbase:
            bases.add(rbase)

        # Obtain class and getter for managed objects
        checks = []
        for uu in c_func.uses:
            obj = self._shared_object(uu.ret_type)
            if obj:
                bases.add(obj)
            else:
                checks.append(uu.name)

        # Ensure that all error codes are set correctly
        error = [-1, "ERR"]
        cxx_rtype = cxx_func.ret_type
        if cxx_rtype.endswith(("int", "size_t")):
            error = ["ERR", "ERR"]
        elif cxx_rtype.endswith("double"):
            error = ["DERR", "DERR"]
        elif "shared_ptr" in cxx_rtype:
            error = ["-2", "ERR"]

        ret = {
            "base": base, "handle": handle,
            "before": before, "buffer": buffer, "after": after,
            "checks": checks, "error": error, "cxx_rbase": rbase,
            "cxx_base": cxx_func.base, "cxx_name": cxx_func.name, "cxx_args": args,
            "cxx_wraps": cxx_wraps,
            "c_func": c_func.name, "c_args": [arg.name for arg in c_func.arglist],
        }
        return ret, bases

    def _scaffold_body(self, c_func: Func, recipe: Recipe) -> tuple[str, set[str]]:
        """Scaffold body of generic CLib function via Jinja."""
        loader = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True,
                             line_comment_prefix="##")
        args, bases = self._crosswalk(c_func, recipe.base)
        args["what"] = recipe.what

        if recipe.code:
            # override auto-generated code
            if "reserved" in recipe.what:
                template = loader.from_string(c_func.wraps)
                after = template.render(cabinets=[kk for kk in self._clib_bases if kk])
            else:
                after = c_func.wraps
            template = loader.from_string(self._templates["clib-custom-code"])
            args["after"] = after.split("\n")

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
            ret_type = c_func.wraps.ret_type
            if "size_t" in ret_type:
                template = loader.from_string(self._templates["clib-size-getter"])
            else:
                template = loader.from_string(self._templates["clib-method"])

        elif recipe.what == "getter":
            ret_type = c_func.wraps.ret_type
            if "void" in ret_type or "vector" in ret_type:
                template = loader.from_string(self._templates["clib-array-getter"])
            elif "size_t" in ret_type:
                template = loader.from_string(self._templates["clib-size-getter"])
            else:
                template = loader.from_string(self._templates["clib-method"])

        elif recipe.what == "setter":
            if "*" in c_func.wraps.arglist[0].p_type:
                template = loader.from_string(self._templates["clib-array-setter"])
            else:
                template = loader.from_string(self._templates["clib-method"])

        else:
            msg = f"{recipe.what!r} not implemented: {c_func.name!r}."
            _LOGGER.critical(msg)
            exit(1)

        macro_module = loader.from_string(self._templates["clib-macros"]).module

        body = template.render(
            render_lines=macro_module.render_lines, get_obj=macro_module.get_obj,
            size_check=macro_module.size_check, **args)
        # remove blank lines left by line comments
        # see https://github.com/pallets/jinja/issues/204
        body = "\n".join(line for line in body.split("\n") if line.strip())
        return body, bases

    def _write_header(self, headers: HeaderFile) -> None:
        """Parse header specification and generate header file."""
        loader = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True)

        filename = headers.output_name(suffix=".h")
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

        preamble = self._config.preambles.get("default")
        preamble_add = self._config.preambles.get(headers.base)
        if preamble_add:
            preamble += "\n" + preamble_add

        guard = filename.name.upper().replace(".", "_")
        t_file = Path(__file__).parent / "template_header.h.j2"
        template = loader.from_string(t_file.read_text(encoding="utf-8"))
        docstring = [escape_token(dd, "Cantera") for dd in headers.docstring]
        output = template.render(
            name=filename.stem, guard=guard, preamble=preamble, prefix=headers.prefix,
            declarations=declarations, base=headers.base, docstring=docstring)

        out = Path(self._out_dir) / "include" / "cantera_clib" / filename.name
        msg = f"  writing {filename.name!r}"
        _LOGGER.info(msg)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(output + "\n", encoding="utf-8")

    def _write_implementation(self, headers: HeaderFile) -> None:
        """Parse header specification and generate implementation file."""
        loader = Environment(loader=BaseLoader, trim_blocks=True, lstrip_blocks=True)

        filename = headers.output_name(suffix=".cpp")
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
        elif headers.base in other:
            other.remove(headers.base)
        includes = []
        for obj in [headers.base] + list(other):
            includes += self._config.includes[obj]

        t_file = Path(__file__).parent / "template_source.cpp.j2"
        template = loader.from_string(t_file.read_text(encoding="utf-8"))
        output = template.render(
            name=filename.stem, implementations=implementations,
            prefix=headers.prefix, base=headers.base, docstring=headers.docstring,
            includes=includes, other=other, str_utils=str_utils)

        out = Path(self._out_dir) / "src" / filename.name
        msg = f"  writing {filename.name!r}"
        _LOGGER.info(msg)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(output + "\n", encoding="utf-8")

    def generate_source(self, headers_files: list[HeaderFile]) -> None:
        """Generate output."""
        self._clib_bases = [f.base for f in headers_files]

        for headers in headers_files:
            self._write_header(headers)
            self._write_implementation(headers)
