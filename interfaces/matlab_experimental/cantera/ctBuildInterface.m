function ctBuildInterface(ctDir, ctIncludeDir, ctGeneratedClibDir, ctLibDir)
    % Build the Cantera MATLAB C++ interface.
    % ctBuildInterface(ctDir, includeDir, ctLibDir) generates and compiles
    % the MATLAB interface for Cantera from the provided header files and a
    % compiled Cantera library.
    %
    % This function automates the full workflow for creating a usable
    % MATLAB package:
    %   1. Calls ctGenerateLibraryDefinitions to create the initial
    %      `definectMatlab.m` interface definition file from the Cantera C++
    %      headers and libraries.
    %   2. Post-processes the generated file with ctEditLibraryDefinitions to
    %      handle edge cases in array argument shapes and cleanup.
    %   3. Uses MATLAB's CLib generation functions to build the interface.
    %   4. Adds the generated interface folder to the MATLAB path.
    %   5. Copies runtime dependencies and saves the updated path.
    %
    % Input arguments:
    %   ctDir              - Root directory of the Cantera source tree.
    %   ctIncludeDir       - Path to the Cantera include directory.
    %   ctGeneratedClibDir - Path to the generated Cantera clib header files.
    %   ctLibDir           - Path to the compiled Cantera library.

    ctDir = string(ctDir);
    ctIncludeDir = string(ctIncludeDir);
    ctGeneratedClibDir = string(ctGeneratedClibDir);
    ctLibDir = string(ctLibDir);
    headerDir = ctGeneratedClibDir;
    includeDir = [ctGeneratedClibDir, ctIncludeDir];
    outputDir = ctDir + "/interfaces/matlab_experimental/cantera";

    % Display the paths (Optional step, for debugging)
    disp("Paths used in the build process:");
    disp("Cantera Root: " + ctDir);
    disp("Header Directory: " + headerDir);
    disp("Include Directory: " + strjoin(includeDir, ", "));
    disp("Cantera Library: " + ctLibDir);
    disp("Output Folder: " + outputDir);

    ctGenerateLibraryDefinitions(includeDir, headerDir, ctLibDir, outputDir);
    ctEditLibraryDefinitions(outputDir);

    % Build C++ Interface
    libDef = feval("definectMatlab");
    build(libDef);
    addpath(fullfile(libDef.OutputFolder,libDef.InterfaceName));
    libDef.copyRuntimeDependencies(Verbose=true);
    savepath;
end
