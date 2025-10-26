function ctBuildInterface(ctToolboxDir, ctIncludeDir, ctLibDir)
    % Build the Cantera MATLAB interface.
    % ctBuildInterface(ctToolboxDir, includeDir, ctLibDir) generates and compiles
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
    %   ctToolboxDir       - Root directory of the Cantera MATLAB toolbox.
    %   ctIncludeDir       - Path to the Cantera include directory.
    %   ctLibDir           - Path to the compiled Cantera library.
    arguments
        ctToolboxDir (1,1) string {mustBeFolder}
        ctIncludeDir (1,1) string {mustBeFolder}
        ctLibDir (1,1) string {mustBeFolder}
    end

    if ~isfolder(ctIncludeDir + "/cantera_clib")
        error('ctBuildInterface:invalidInclude', ...
            'Invalid include folder. "%s" must have "cantera_lib" subfolder ', ...
            ctIncludeDir);
    end

    ctToolboxDir = string(ctToolboxDir);
    ctIncludeDir = string(ctIncludeDir);
    ctLibDir = string(ctLibDir);
    outputDir = ctToolboxDir + "/cantera";

    % Display the paths (Optional step, for debugging)
    disp("Paths used in the build process:");
    disp("Cantera Toolbox Root: " + ctToolboxDir);
    disp("Include Directory: " + ctIncludeDir);
    disp("Cantera Library: " + ctLibDir);
    disp("Output Folder: " + outputDir);

    ctGenerateLibraryDefinitions(ctIncludeDir, ctLibDir, outputDir);
    ctEditLibraryDefinitions(outputDir);

    % Build C++ Interface
    libDef = feval("definectMatlab");
    build(libDef);
    addpath(fullfile(libDef.OutputFolder,libDef.InterfaceName));
    libDef.copyRuntimeDependencies(Verbose=true);
    delete(fullfile(outputDir, 'definectMatlab.m'));
    delete(fullfile(outputDir, 'ctMatlabData.xml'));
    savepath;
end
