function ctLoad()
    % ctLoad
    % Load the Cantera C Library into Memory

    paths = ctPaths();

    if any(cellfun(@isempty, {paths.libPath, paths.includePath, paths.toolboxPath}))
        error('ctLoad:MissingPath', ...
              'Library, header, and toolbox paths must be configured with ctPaths(libPath,includePath, toolboxPath).');
    end

    if ispc
        libName = 'cantera_shared.dll';
    elseif ismac
        libName = 'libcantera_shared.dylib';
    elseif isunix
        libName = 'libcantera_shared.so';
    else
        error('ctLoad:UnsupportedPlatform', 'Operating system not supported.');
    end

    fullLibPath = fullfile(paths.libPath, libName);
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
