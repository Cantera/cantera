function ctLoad()
    % Load the Cantera C Library into Memory

    if ispc
        ctName = '/bin/cantera_shared.dll';
    elseif ismac
        ctName = '/Lib/libcantera_shared.dylib';
    elseif isunix
        ctName = '/lib/libcantera_shared.so';
    else
        error('Operating System Not Supported!');
        return;
    end

    root = ctRoot;
    if ispc
        root = [ctRoot, '/Library'];
    end

    if ~libisloaded(ctLib)
        [~, warnings] = loadlibrary([root, ctName], ...
                                    [root, '/include/cantera/clib/ctmatlab.h'], ...
                                    'includepath', [root, '/include'], ...
                                    'addheader', 'ct', 'addheader', 'ctfunc', ...
                                    'addheader', 'ctmultiphase', 'addheader', ...
                                    'ctonedim', 'addheader', 'ctreactor', ...
                                    'addheader', 'ctrpath', 'addheader', 'ctsurf');
    end

    disp(sprintf('Cantera %s is ready for use.', ctVersion))

end
