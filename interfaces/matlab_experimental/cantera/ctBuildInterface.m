function ctBuildInterface(mode)
    % Mode determines how the paths are set. By default it's set to 'auto'
    % and the environment variable CANTERA_ROOT is used to determine the
    % other paths. This mode is used when the user compiled Cantera from
    % source, so that all the necessary files are already in their correct
    % locations. If the user obtained the MATLAB toolbox from FileExchange
    % and the header files and compiled library files are downloaded by the
    % MLTBX setup script, the 'FileExchange' mode is used.

    if nargin < 1
        mode = 'auto';
    end

    if strcmp(mode, 'auto')
        canteraRoot = string(getenv('CANTERA_ROOT'));
        includeDir = canteraRoot + "/include";
        headerDir = includeDir + "/cantera/clib";
        ctLibDir = canteraRoot + "/build/lib";
        outputDir = canteraRoot + "/interfaces/matlab_experimental/cantera";
    elseif strcmp(mode, 'FileExchange')
        canteraRoot = string(pwd);
        includeDir = canteraRoot + "/include";
        headerDir = includeDir;
        ctLibDir = canteraRoot + "/lib";
        outputDir = canteraRoot + "/matlab_experimental/cantera";
    else
        error('Invalid mode. Use "auto" or "FileExchange".');
    end

    % Display the paths (Optional step, for debugging)
    disp("Paths used in the build process:");
    disp('Cantera Root: ' + canteraRoot);
    disp('Include Directory: ' + includeDir);
    disp('Header Directory: '+ headerDir);
    disp('Cantera Library: ' + ctLibDir);
    disp('Output Folder: ' + outputDir);

    ctGenerateLibraryDefinitions(includeDir, headerDir, ctLibDir, outputDir);
    ctEditLibraryDefinitions(outputDir);

    % Print out the library definitions
    libDef = feval("definectMatlab");

    build(libDef);
    addpath(fullfile(libDef.OutputFolder,libDef.InterfaceName));
    libDef.copyRuntimeDependencies(Verbose=true);
    savepath;
end
