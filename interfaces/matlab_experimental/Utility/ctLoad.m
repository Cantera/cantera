function ctLoad(paths)
    % ctLoad
    % Load the Cantera C Library into Memory
    %
    %   ctLoad()          uses paths from ctPaths(), default
    %   ctLoad(paths)     uses explicitly provided paths (only used by unit tests)

    if nargin < 1
        paths = ctPaths();
    end

    if any(cellfun(@isempty, {paths.libPath, paths.includePath, paths.toolboxPath}))
        error('ctLoad:MissingPath', ...
              'Library, header, and toolbox paths must be configured with `ctPaths`');
    end

    if ispc
        pattern = fullfile(paths.libPath, [ctLib, '*.dll']);
    elseif ismac
        pattern = fullfile(paths.libPath, [ctLib, '*.dylib']);
    elseif isunix
        pattern = fullfile(paths.libPath, [ctLib, '*.so']);
    else
        error('ctLoad:UnsupportedPlatform', 'Operating system not supported.');
    end

    libFiles = dir(pattern);
    if isempty(libFiles)
        error('ctLoad:LibraryNotFound', ...
              'No matching Cantera library found in %s (pattern: %s)', paths.libPath, pattern);
    end

    libFile = libFiles(1).name;
    fullLibPath = fullfile(paths.libPath, libFile);
    fullHeaderPath = fullfile(paths.includePath, 'cantera', 'clib', 'ctmatlab.h');

    if ~libisloaded(ctLib)
        [~, warnings] = loadlibrary(fullLibPath, fullHeaderPath, ...
                                    'includepath', paths.includePath, ...
                                    'addheader', 'ct', ...
                                    'addheader', 'ctfunc', ...
                                    'addheader', 'ctmultiphase', ...
                                    'addheader', 'ctonedim', ...
                                    'addheader', 'ctreactor', ...
                                    'addheader', 'ctrpath', ...
                                    'addheader', 'ctsurf');
    end

    fprintf('Cantera %s is ready for use.\n', ctVersion);

end
