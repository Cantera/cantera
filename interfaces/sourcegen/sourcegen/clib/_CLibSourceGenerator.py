"""Generator for CLib source files."""

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
import logging
from dataclasses import dataclass

from jinja2 import Environment, BaseLoader

from .._dataclasses import HeaderFile, Param, ArgList, CFunc, Recipe
from .._SourceGenerator import SourceGenerator
from .._helpers import with_unpack_iter


_LOGGER = logging.getLogger()


@dataclass(frozen=True)
@with_unpack_iter
class Config:
    """Provides configuration info for the CLibSourceGenerator class"""

    ret_type_crosswalk: dict[str, str]  #: Return type cross-walks

    prop_type_crosswalk: dict[str, str]  #: Parameter type cross-walks

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
        elif isinstance(c_func.implements, Param):
            implements += f": {c_func.implements.long_str()}"
        elif isinstance(c_func.implements, str):
            implements += ": custom code"
        block = template.render(
            brief=c_func.brief,
            params=[param(par) for par in c_func.arglist],
            returns=c_func.returns, implements=implements,
            relates=[f"{uu.base}::{uu.name}()" for uu in c_func.uses])
        return self._javadoc_comment(block)

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
        elif isinstance(cxx_member, str):
            cxx_func = CFunc("void", "dummy", ArgList([]), "", None, "", "base")
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
                elif "double" in cxx_type or "void" in cxx_type:
                    # can pass values directly
                    buffer = []
                elif "int" in cxx_type or "size_t" in cxx_type:
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

        cxx_implements = ""
        if cxx_func.name != "dummy":
            cxx_implements = cxx_func.short_declaration()

        ret = {
            "base": base, "handle": handle, "lines": lines, "buffer": buffer,
            "shared": shared, "checks": checks, "error": error, "cxx_rbase": cxx_rbase,
            "cxx_base": cxx_func.base, "cxx_name": cxx_func.name, "cxx_args": args,
            "cxx_implements": cxx_implements,
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
            args["lines"] = c_func.implements.split("\n")

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
        t_file = Path(__file__).parent / "template_header.h.in"
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

        t_file = Path(__file__).parent / "template_source.cpp.in"
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

    def generate_source(self, headers_files: list[HeaderFile]) -> None:
        """Generate output."""
        self._clib_bases = [f.base for f in headers_files]

        for headers in headers_files:
            self._write_header(headers)
            self._write_implementation(headers)
