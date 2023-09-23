clear all

% Copy library to test folder
ctTestPath;

if ispc
    ctName = '/build/lib/cantera_shared.dll';
elseif ismac
    ctName = '/build/lib/libcantera_shared.dylib';
elseif isunix
    ctName = '/build/lib/libcantera_shared.so';
end
% Load Cantera
if ~libisloaded('libcantera_shared')
    [~, warnings] = loadlibrary([cantera_root, ctName], ...
                                [cantera_root, '/include/cantera/clib/ctmatlab.h'], ...
                                'includepath', [cantera_root, '/include'], ...
                                'addheader', 'ct', 'addheader', 'ctfunc', ...
                                'addheader', 'ctmultiphase', 'addheader', ...
                                'ctonedim', 'addheader', 'ctreactor', ...
                                'addheader', 'ctrpath', 'addheader', 'ctsurf');
end

disp('Cantera is loaded for test');
