clear all

% Copy library to test folder
ctTestPath;

% Load Cantera
rootDir = getenv('CANTERA_ROOT');

if ispc
    ctName = '/test/matlab_experimental/cantera_shared.dll';
elseif ismac
    ctname = '/test/matlab_experimental/libcantera_shared.dylib';
elseif isunix
    ctname = '/test/matlab_experimental/libcantera_shared.so';
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
