function LoadCantera
    addpath('Class', 'Utility', 'PresetMixtures', 'PresetReactors', ...
            'PresetObjects', '1D');
    if ispc
        ctname = 'cantera_shared.dll';
    elseif ismac
        ctname = 'libcantera_shared.dylib';
    elseif isunix
        ctname = 'libcantera_shared.so';

    else
        error('Operating System Not Supported!');
        return;
    end
    if ~libisloaded(ct)
        load('Utility/cantera_root.mat');
        [~,warnings] = loadlibrary([cantera_root '/lib/' ctname], ...
                                   [cantera_root '/include/cantera/clib/ctmatlab.h'], ...
                                   'includepath', [cantera_root '/include'], ...
                                   'addheader','ct','addheader','ctfunc', ...
                                   'addheader','ctmultiphase','addheader', ...
                                   'ctonedim','addheader','ctreactor', ...
                                   'addheader','ctrpath','addheader','ctsurf', ...
                                   'addheader','ctxml');
        example_dir = [cantera_root, '/samples/matlab_new'];
        addpath(example_dir);
    end
    disp('Cantera is ready for use');
end
