% Load Cantera without changing rootDir
rootDir = '/home/ssun30/anaconda3/envs/ct-matlab';
ctName = '/lib/libcantera_shared.so';
% Load Cantera
if ~libisloaded('libcantera_shared')
    [~, warnings] = loadlibrary([rootDir, ctName], ...
                                [rootDir, '/include/cantera/clib/ctmatlab.h'], ...
                                'includepath', [rootDir '/include'], ...
                                'addheader', 'ct', 'addheader', 'ctfunc', ...
                                'addheader', 'ctmultiphase', 'addheader', ...
                                'ctonedim', 'addheader', 'ctreactor', ...
                                'addheader', 'ctrpath', 'addheader', 'ctsurf');
end

disp('Cantera is loaded for test');

% Run all tests
runtests('ctMatlabTestThermo');

% Unload Cantera
unloadlibrary('libcantera_shared');
disp('Cantera has been unloaded');