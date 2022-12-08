% Load the Cantera C Library into the Memory
%
if ispc
    ctname = '/Lib/cantera_shared.dll';
elseif ismac
    ctname = '/Lib/libcantera_shared.dylib';
elseif isunix
    ctname = '/lib/libcantera_shared.so';
else
    error('Operating System Not Supported!');
    return;
end

if ~libisloaded(ct)
    [~, warnings] = loadlibrary([cantera_root, ctname], ...
                                [cantera_root, '/include/cantera/clib/ctmatlab.h'], ...
                                'includepath', [cantera_root '/include'], ...
                                'addheader', 'ct', 'addheader', 'ctfunc', ...
                                'addheader', 'ctmultiphase', 'addheader', ...
                                'ctonedim', 'addheader', 'ctreactor', ...
                                'addheader', 'ctrpath', 'addheader', 'ctsurf');
end

ct_ver = canteraVersion;
sprintf('%s is ready for use.', ct_ver);
clear all
