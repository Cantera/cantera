#!/usr/bin/env julia
# generate_bindings.jl
#
# Reproducible generator for the low-level Julia bindings to Cantera's
# generated CLib API (the `cantera_clib/ct*.h` headers).
#
# This script parses the CLib C headers and emits one Julia file per header
# under `src/generated/`, containing a thin `ccall` wrapper for every exported
# function.  The generated wrappers are INTERNAL: the public Julia
# API lives in the hand-written `src/*.jl` files and must never be generated.
#
# Design goals:
#   * deterministic output (headers processed in sorted order, functions in
#     declaration order);
#   * no hand-editing of generated files (regenerate instead);
#   * minimal, mechanical C->Julia type mapping so that churn in the
#     experimental CLib only requires re-running this script.
#
# Usage:
#   julia generate/generate_bindings.jl [--headers <dir>] [--out <dir>]
#
# If `--headers` is omitted the script looks for the headers in, in order:
#   1. $CANTERA_CLIB_INCLUDE
#   2. a sibling Cantera build tree (../clib/include/cantera_clib)
#   3. an active conda environment ($CONDA_PREFIX/include/cantera_clib)
#
# The headers follow a very regular contract (see interfaces/clib):
#   * every function returns `int32_t` or `double`;
#   * `int32_t` return is a handle, count, index, status (0/-1) or, for string
#     getters, the required buffer length; `-1`/`ERR` signals an exception;
#   * `double` return is a scalar property; `DERR` signals an exception;
#   * strings in:  `const char*`;  strings out: `int32_t bufLen, char* buf`;
#   * arrays in:   `int32_t nLen, const double* n`;
#   * arrays out:  `int32_t bufLen, double* buf`.

const CTYPE_MAP = Dict(
    "int32_t"        => "Int32",
    "const int32_t"  => "Int32",
    "int64_t"        => "Int64",
    "const int64_t"  => "Int64",
    "double"         => "Float64",
    "const double"   => "Float64",
    "const char*"    => "Cstring",
    "char*"          => "Ptr{UInt8}",
    "const double*"  => "Ptr{Float64}",
    "double*"        => "Ptr{Float64}",
    "const int32_t*" => "Ptr{Int32}",
    "int32_t*"       => "Ptr{Int32}",
    "LogCallback"    => "Ptr{Cvoid}",
)

"Strip C `/* */` and `//` comments from a source string."
function strip_comments(src::AbstractString)
    src = replace(src, r"/\*.*?\*/"s => " ")
    src = replace(src, r"//[^\n]*" => " ")
    return src
end

"Normalize a C type fragment (collapse whitespace, keep `const` and `*`)."
function normalize_ctype(s::AbstractString)
    s = strip(s)
    s = replace(s, r"\s*\*" => "*")   # `char *` -> `char*`
    s = replace(s, r"\s+" => " ")
    return s
end

struct CArg
    ctype::String
    name::String
end

struct CFunc
    ret::String
    name::String
    args::Vector{CArg}
end

"Split a comma-separated C argument list at top level (no nested parens here)."
function split_args(s::AbstractString)
    s = strip(s)
    (isempty(s) || s == "void") && return String[]
    return strip.(split(s, ','))
end

"Parse a single argument fragment like `const double* x` into (ctype, name)."
function parse_arg(frag::AbstractString)
    frag = normalize_ctype(frag)
    m = match(r"^(.*?)(\b[A-Za-z_]\w*)$", frag)  # last identifier is the name
    m === nothing && error("cannot parse argument: `$frag`")
    ctype = normalize_ctype(m.captures[1])
    name = m.captures[2]
    return CArg(ctype, name)
end

"Parse all function declarations from one (comment-stripped) header body."
function parse_header(src::AbstractString)
    src = strip_comments(src)
    funcs = CFunc[]
    for m in eachmatch(r"\b(int32_t|int64_t|double|void)\s+([A-Za-z_]\w*)\s*\(([^;{]*)\)\s*;"s, src)
        ret, name, arglist = m.captures
        args = CArg[parse_arg(a) for a in split_args(arglist)]
        push!(funcs, CFunc(ret, name, args))
    end
    return funcs
end

function julia_type(ctype::AbstractString)
    haskey(CTYPE_MAP, ctype) && return CTYPE_MAP[ctype]
    error("unmapped C type: `$ctype` (extend CTYPE_MAP in generate_bindings.jl)")
end

function emit_func(io, f::CFunc)
    jlret = julia_type(f.ret)
    argnames = [a.name for a in f.args]
    argtypes = [julia_type(a.ctype) for a in f.args]
    sig = join(argnames, ", ")
    tuple = isempty(argtypes) ? "()" : "(" * join(argtypes, ", ") * ",)"
    call_args = isempty(argnames) ? "" : ", " * join(argnames, ", ")
    println(io, "function $(f.name)($sig)")
    println(io, "    ccall((:$(f.name), libcantera[]), $jlret, $tuple$call_args)")
    println(io, "end")
end

const HEADER_ORDER = ["ct", "ctsol", "ctthermo", "ctkin", "cttrans",
                      "ctrxn", "ctreactor", "ctreactornet", "ctonedim",
                      "ctdomain", "ctmix", "ctfunc", "ctrdiag", "ctconnector"]

function find_headers()
    if haskey(ENV, "CANTERA_CLIB_INCLUDE")
        return ENV["CANTERA_CLIB_INCLUDE"]
    end
    candidates = String[
        joinpath(@__DIR__, "..", "..", "clib", "include", "cantera_clib"),
    ]
    haskey(ENV, "CONDA_PREFIX") &&
        push!(candidates, joinpath(ENV["CONDA_PREFIX"], "include", "cantera_clib"))
    for c in candidates
        isdir(c) && return c
    end
    error("could not locate cantera_clib headers; set CANTERA_CLIB_INCLUDE")
end

function main(args)
    headers_dir = nothing
    out_dir = joinpath(@__DIR__, "..", "src", "generated")
    i = 1
    while i <= length(args)
        if args[i] == "--headers"
            headers_dir = args[i+1]; i += 2
        elseif args[i] == "--out"
            out_dir = args[i+1]; i += 2
        else
            error("unknown argument: $(args[i])")
        end
    end
    headers_dir = headers_dir === nothing ? find_headers() : headers_dir
    mkpath(out_dir)

    available = filter(h -> isfile(joinpath(headers_dir, "$h.h")), HEADER_ORDER)
    isempty(available) && error("no known headers found in $headers_dir")

    generated_files = String[]
    total = 0
    for h in available
        src = read(joinpath(headers_dir, "$h.h"), String)
        funcs = parse_header(src)
        outfile = joinpath(out_dir, "lib$h.jl")
        open(outfile, "w") do io
            println(io, "# This file was generated by generate/generate_bindings.jl.")
            println(io, "# Source header: cantera_clib/$h.h")
            println(io, "# DO NOT EDIT: regenerate with `julia generate/generate_bindings.jl`.")
            println(io)
            for (i, f) in enumerate(funcs)
                i > 1 && println(io)
                emit_func(io, f)
            end
        end
        push!(generated_files, "lib$h.jl")
        total += length(funcs)
        println("  lib$h.jl: $(length(funcs)) functions")
    end

    # Manifest that LibCantera.jl includes.
    open(joinpath(out_dir, "_manifest.jl"), "w") do io
        println(io, "# Generated manifest of low-level binding files (see generate_bindings.jl).")
        println(io, "const GENERATED_FILES = [")
        for f in generated_files
            println(io, "    \"$f\",")
        end
        println(io, "]")
    end
    println("Generated $total wrappers across $(length(available)) headers into $out_dir")
end

main(ARGS)
