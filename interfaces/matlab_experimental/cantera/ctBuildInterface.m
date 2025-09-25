function ctBuildInterface(ctDir, includeDir, ctLibDir)
    % Builds the Cantera MATLAB C++ Interface from the Cantera CLib
    % header files and the compiled Cantera library.

    ctDir = string(ctDir);
    includeDir = string(includeDir);
    ctLibDir = string(ctLibDir);
    headerDir = includeDir + "/cantera/clib";
    outputDir = ctDir + "/interfaces/matlab_experimental/cantera";

    % Display the paths (Optional step, for debugging)
    disp("Paths used in the build process:");
    disp("Cantera Root: " + ctDir);
    disp("Include Directory: " + includeDir);
    disp("Header Directory: "+ headerDir);
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
