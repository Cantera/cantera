function ctGenerateLibraryDefinitions(includeDir, ctLibDir, outputDir)
    % Generate the MATLAB interface definition file.
    % ctGenerateLibraryDefinitions(includeDir, ctLibDir, outputDir)
    % creates the `definectMatlab.m` interface definition file for the
    % Cantera C library. This file is later compiled into the MATLAB interface.
    %
    %   The function:
    %     1. Collects all header files under headerDir.
    %     2. Locates the compiled Cantera shared library in ctLibDir.
    %     3. Configures the MATLAB C++ compiler (MEX).
    %     4. Calls clibgen to produce the `definectMatlab.m` file in outputDir.

    headerFiles = dir(fullfile(includeDir, 'cantera_clib', '*.h'));
    headerPaths = string({headerFiles.name});
    headerPaths = fullfile({headerFiles.folder}, headerPaths);

    % Get path for the shared library file
    libraries = ctLib(ctLibDir);
    disp("Using shared library: " + libraries);

    if isMATLABReleaseOlderThan("R2024a")
        nameArg = "PackageName";
    else
        nameArg = "InterfaceName";
    end

    overwriteExistingDefinitionFiles = true;

    % Set up C++ compiler
    mex -setup cpp

    % Generate definition file for C++ library
    clibgen.generateLibraryDefinition(headerPaths, ...
        "IncludePath", includeDir, ...
        "Libraries", libraries, ...
        "OutputFolder", outputDir, ...
        nameArg, "ctMatlab", ...
        "OverwriteExistingDefinitionFiles", overwriteExistingDefinitionFiles, ...
        "CLinkage", true, ...
        "TreatObjectPointerAsScalar", true, ...
        "TreatConstCharPointerAsCString", true, ...
        "ReturnCArrays", false, ...
        "Verbose", true);
end
