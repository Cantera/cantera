# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Library-name matching only: needs no built Cantera library, so it can also be run
# on its own with `julia test/test_library_names.jl`.

using Test

module LibraryNames
include(joinpath(@__DIR__, "..", "src", "library_names.jl"))
end

const LN = LibraryNames

@testset "Shared library names" begin
    @testset "Windows (gh-2192)" begin
        # Files in build\lib after an MSVC build
        msvc = ["cantera.lib", "cantera_shared.dll", "cantera_shared.lib",
                "cantera_shared.exp", "cantera_python3_15.dll"]
        @test LN.pick_library(msvc; windows=true) == "cantera_shared.dll"
        @test LN.pick_library(["cantera.lib", "cantera.dll"]; windows=true) == "cantera.dll"
        # MinGW keeps the lib prefix
        @test LN.pick_library(["libcantera_shared.dll", "libcantera.a"];
                              windows=true) == "libcantera_shared.dll"
        # Case-insensitive file system
        @test LN.pick_library(["Cantera_Shared.DLL"]; windows=true) == "Cantera_Shared.DLL"
        # Never the static/import library or a Python extension module
        @test LN.pick_library(["cantera.lib", "cantera_python3_15.dll"];
                              windows=true) === nothing
        @test LN.pick_library(["cantera_shared.lib"]; windows=true) === nothing
        @test LN.fallback_library_name(; windows=true) == "cantera_shared"
        dirs = LN.library_search_dirs(; windows=true,
                                      env=Dict("CONDA_PREFIX" => "C:\\envs\\ct"))
        @test dirs[1] == joinpath("C:\\envs\\ct", "Library", "bin")
        @test !("/usr/lib" in dirs)
    end

    @testset "Linux and macOS (unchanged)" begin
        @test LN.pick_library(["libcantera_shared.so", "libcantera.a"];
                              windows=false) == "libcantera_shared.so"
        @test LN.pick_library(["libcantera_shared.so.3.2.0"];
                              windows=false) == "libcantera_shared.so.3.2.0"
        @test LN.pick_library(["libcantera_shared.3.dylib"];
                              windows=false) == "libcantera_shared.3.dylib"
        @test LN.pick_library(["cantera_shared.dll"]; windows=false) === nothing
        @test LN.fallback_library_name(; windows=false) == "libcantera"
        @test "/usr/lib" in LN.library_search_dirs(; windows=false, env=Dict{String,String}())
    end
end
