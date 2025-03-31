function ctGenerateLibraryDefinitions(includeDir, headerDir, ctLibDir, outputDir)
    % Generate paths for header files
    headerFiles = dir(fullfile(headerDir, '**', '*.h'));
    headerPaths = string({headerFiles.name});
    headerPaths = fullfile({headerFiles.folder}, headerPaths);

    % Get path for the shared library file
    if ispc
        libraries = ctLibDir + "/cantera_shared.dll";
    elseif islinux
        libraries = ctLibDir + "/libcantera_shared.so";
    elseif ismac
        libraries = ctLibDir + "/libcantera_shared.dylib";
    else
        error('Operating System Not Supported');
    end

    overwriteExistingDefinitionFiles = true;

    % Set up C++ compiler
    mex -setup cpp

    % Generate definition file for C++ library
    clibgen.generateLibraryDefinition(headerPaths, ...
        "IncludePath",includeDir, ...
        "Libraries",libraries, ...
        "OutputFolder",outputDir, ...
        "InterfaceName","ctMatlab", ...
        "OverwriteExistingDefinitionFiles",overwriteExistingDefinitionFiles, ...
        "CLinkage",true, ...
        "TreatObjectPointerAsScalar",true, ...
        "TreatConstCharPointerAsCString",true, ...
        "ReturnCArrays",false, ...
        "Verbose",true);
end
