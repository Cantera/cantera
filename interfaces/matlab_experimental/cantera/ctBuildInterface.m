function ctBuildInterface()
    % Builds the Cantera MATLAB C++ Interface from the Cantera CLib
    % header files and the compiled Cantera library.
    % The environment variable CANTERA_ROOT is used to obtain the paths.

    canteraRoot = string(getenv('CANTERA_ROOT'));
    includeDir = canteraRoot + "/include";
    headerDir = includeDir + "/cantera/clib";
    ctLibDir = canteraRoot + "/build/lib";
    outputDir = canteraRoot + "/interfaces/matlab_experimental/cantera";

    % Display the paths (Optional step, for debugging)
    disp("Paths used in the build process:");
    disp('Cantera Root: ' + canteraRoot);
    disp('Include Directory: ' + includeDir);
    disp('Header Directory: '+ headerDir);
    disp('Cantera Library: ' + ctLibDir);
    disp('Output Folder: ' + outputDir);

    ctGenerateLibraryDefinitions(includeDir, headerDir, ctLibDir, outputDir);
    ctEditLibraryDefinitions(outputDir);

    % Build C++ Interface
    libDef = feval("definectMatlab");
    build(libDef);
    addpath(fullfile(libDef.OutputFolder,libDef.InterfaceName));
    libDef.copyRuntimeDependencies(Verbose=true);
    savepath;
end
