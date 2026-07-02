"""
    LibCantera

Internal module holding the low-level `ccall` bindings to Cantera's generated
CLib API.  Everything here is an implementation detail: user code should use the
API exported from the top-level [`Cantera`](@ref) module instead of the
raw handle-based functions defined here.

The wrapper functions in `generated/` are produced by
`generate/generate_bindings.jl` from the `cantera_clib/*.h` headers and must not
be hand-edited.  This module is responsible only for:

  * locating the Cantera shared library (`libcantera`);
  * locating Cantera's data directory (for `gri30.yaml` and friends);
  * `include`-ing the generated binding files.
"""
module LibCantera

export libcantera

"""
    find_library() -> String

Locate the Cantera shared library.  Search order:

 1. `\$CANTERA_LIBRARY_PATH` (a directory or a full path to the library);
 2. `\$CONDA_PREFIX/lib` for a conda-installed Cantera;
 3. a sibling Cantera build tree (`../../../build/lib`);
 4. the system loader's default search path (`libcantera`).

Throws an informative error if nothing is found.
"""
function find_library()
    # Installed names first, then the in-tree build-tree name (libcantera_shared).
    names = ["libcantera.so", "libcantera.dylib", "libcantera.dll",
             "libcantera.so.3", "libcantera_shared.so", "libcantera_shared.dylib",
             "libcantera_shared.dll", "cantera"]

    function probe(dir)
        isdir(dir) || return nothing
        for n in names
            p = joinpath(dir, n)
            isfile(p) && return p
        end
        # also accept versioned sonames like libcantera.so.3.2.0
        for f in readdir(dir)
            startswith(f, "libcantera") && occursin(".so", f) && return joinpath(dir, f)
            startswith(f, "libcantera") && endswith(f, ".dylib") &&
                return joinpath(dir, f)
        end
        return nothing
    end

    if haskey(ENV, "CANTERA_LIBRARY_PATH")
        p = ENV["CANTERA_LIBRARY_PATH"]
        isfile(p) && return p
        hit = probe(p)
        hit === nothing || return hit
    end

    dirs = String[]
    haskey(ENV, "CONDA_PREFIX") && push!(dirs, joinpath(ENV["CONDA_PREFIX"], "lib"))
    push!(dirs, normpath(joinpath(@__DIR__, "..", "..", "..", "build", "lib")))
    append!(dirs, ["/usr/local/lib", "/usr/lib"])
    for d in dirs
        hit = probe(d)
        hit === nothing || return hit
    end

    # Fall back to the loader's search path; dlopen will error if truly missing.
    return "libcantera"
end

# Resolved at runtime in `__init__` (below), NOT baked into the precompile image
# — otherwise the discovered path would be frozen at precompile time and later
# changes to CANTERA_LIBRARY_PATH would be silently ignored.  Generated `ccall`
# wrappers dereference this Ref as `libcantera[]`.
const libcantera = Ref{String}("libcantera")

function __init__()
    libcantera[] = find_library()
    return nothing
end

# ---- generated low-level bindings -------------------------------------------
include("generated/_manifest.jl")
for f in GENERATED_FILES
    include(joinpath("generated", f))
end

"""
    default_data_directory() -> Union{String,Nothing}

Best-effort discovery of Cantera's bundled data directory (containing
`gri30.yaml`).  Used to register the directory with Cantera at load time so that
mechanism files can be found by name.
"""
function default_data_directory()
    dirs = String[]
    haskey(ENV, "CANTERA_DATA") && push!(dirs, ENV["CANTERA_DATA"])
    if haskey(ENV, "CONDA_PREFIX")
        push!(dirs, joinpath(ENV["CONDA_PREFIX"], "share", "cantera", "data"))
    end
    push!(dirs, normpath(joinpath(@__DIR__, "..", "..", "..", "data")))
    for d in dirs
        isdir(d) && isfile(joinpath(d, "gri30.yaml")) && return d
    end
    return nothing
end

end # module LibCantera
