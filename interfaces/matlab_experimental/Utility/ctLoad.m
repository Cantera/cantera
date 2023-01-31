function ctLoad()
    % Load the Cantera C Library into Memory
    %
    if ispc
        ctName = '/Lib/cantera_shared.dll';
    elseif ismac
        ctName = '/Lib/libcantera_shared.dylib';
    elseif isunix
        ctName = '/lib/libcantera_shared.so';
    else
        error('Operating System Not Supported!');
        return;
    end

    if ~libisloaded(ctSharedLibrary)
        [~, warnings] = loadlibrary([ctRoot, ctName], ...
                                    [ctRoot, '/include/cantera/clib/ctmatlab.h'], ...
                                    'includepath', [ctRoot '/include'], ...
                                    'addheader', 'ct', 'addheader', 'ctfunc', ...
                                    'addheader', 'ctmultiphase', 'addheader', ...
                                    'ctonedim', 'addheader', 'ctreactor', ...
                                    'addheader', 'ctrpath', 'addheader', 'ctsurf');
    end

    disp(sprintf('Cantera %s is ready for use.', ctVersion))

end
