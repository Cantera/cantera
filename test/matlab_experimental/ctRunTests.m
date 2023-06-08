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

% Run all tests
results1 = runtests('ctTestThermo')

% Unload Cantera and remove temporary library file
unloadlibrary('libcantera_shared');
delete([rootDir, ctName]);

disp('Cantera has been unloaded');