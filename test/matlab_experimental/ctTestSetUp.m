clear all

% Copy library to test folder
ctTestPath;

% Load Cantera
rootDir = getenv('CANTERA_ROOT');

if ispc
    ctName = '/build/lib/cantera_shared.dll';
elseif ismac
    ctName = '/build/lib/libcantera_shared.dylib';
elseif isunix
    ctName = '/build/lib/libcantera_shared.so';
end
% Load Cantera
if ~libisloaded('libcantera_shared')
    [~, warnings] = loadlibrary([rootDir, ctName], ...
                                [rootDir, '/include/cantera/clib/ctmatlab.h'], ...
                                'includepath', [rootDir, '/include'], ...
                                'addheader', 'ct', 'addheader', 'ctfunc', ...
                                'addheader', 'ctmultiphase', 'addheader', ...
                                'ctonedim', 'addheader', 'ctreactor', ...
                                'addheader', 'ctrpath', 'addheader', 'ctsurf');
end

disp('Cantera is loaded for test');
