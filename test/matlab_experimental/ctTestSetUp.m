clear all

% Copy library to test folder
ctTestPath;

% Load Cantera
rootDir = fullfile(pwd);
ctName = '/test/matlab_experimental/libcantera_shared.so';
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
