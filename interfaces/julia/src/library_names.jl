# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Platform-specific names of the Cantera shared library. Kept free of any `ccall`
# so that it can be tested without a built library.

"""
    library_names(; windows=Sys.iswindows()) -> Vector{String}

File names of the Cantera shared library, in the order they are tried.

On Windows, DLLs have no `lib` prefix: an MSVC build produces
`cantera_shared.dll` (or `cantera.dll` with `renamed_shared_libraries=n`),
installed to `bin` or, for conda, `Library\\bin`. MinGW builds use a `lib`
prefix. Elsewhere the library is `libcantera[_shared].so` or `.dylib`.
"""
function library_names(; windows::Bool = Sys.iswindows())
    if windows
        return ["cantera_shared.dll", "cantera.dll",
                "libcantera_shared.dll", "libcantera.dll"]
    end
    return ["libcantera.so", "libcantera.dylib", "libcantera.so.3",
            "libcantera_shared.so", "libcantera_shared.dylib", "cantera"]
end

"""
    pick_library(files; windows=Sys.iswindows()) -> Union{String,Nothing}

Choose the Cantera shared library from the file names `files` found in one
directory, or return `nothing`. Exact names from [`library_names`](@ref) win;
otherwise versioned names are accepted (`libcantera.so.3.2.0`,
`libcantera_shared.3.dylib`, `cantera_shared-3.dll`). On Windows, import or
static libraries (`.lib`) and Python extension modules
(`cantera_python3_15.dll`) are never chosen.
"""
function pick_library(files::AbstractVector{<:AbstractString};
                      windows::Bool = Sys.iswindows())
    # Windows file names are case-insensitive
    canon(f) = windows ? lowercase(f) : f
    present = Dict(canon(f) => String(f) for f in files)
    for n in library_names(; windows)
        haskey(present, n) && return present[n]
    end
    for f in files
        if windows
            occursin(r"^(lib)?cantera(_shared)?[-.]?[0-9.]*\.dll$", lowercase(f)) &&
                return String(f)
        else
            startswith(f, "libcantera") && occursin(".so", f) && return String(f)
            startswith(f, "libcantera") && endswith(f, ".dylib") && return String(f)
        end
    end
    return nothing
end

"""
    library_search_dirs(; windows=Sys.iswindows(), env=ENV) -> Vector{String}

Directories searched (after `\$CANTERA_LIBRARY_PATH`) for the shared library.
"""
function library_search_dirs(; windows::Bool = Sys.iswindows(), env = ENV)
    dirs = String[]
    if haskey(env, "CONDA_PREFIX")
        # conda puts DLLs in Library\bin on Windows
        windows && push!(dirs, joinpath(env["CONDA_PREFIX"], "Library", "bin"))
        push!(dirs, joinpath(env["CONDA_PREFIX"], "lib"))
    end
    push!(dirs, normpath(joinpath(@__DIR__, "..", "..", "..", "build", "lib")))
    windows || append!(dirs, ["/usr/local/lib", "/usr/lib"])
    return dirs
end

"""Name passed to the system loader when no file was found."""
fallback_library_name(; windows::Bool = Sys.iswindows()) =
    windows ? "cantera_shared" : "libcantera"
