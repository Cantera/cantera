% testpath.m
% Set up environment for testing the Cantera Matlab interface
% from within the Cantera source tree. Run this file from the
% root of the Cantera source tree, for example:
%
%    cd ~/src/cantera
%    run interfaces/matlab/testpath.m

% Unload the mex file so copying the DLL will work
clear ctmethods

% get the list of directories on the Matlab path
dirs = regexp(path, ['([^' pathsep ']*)'], 'match');

% if 'cantera' is already in the path, remove it
for i = 1:length(dirs)
    if strfind(dirs{i}, 'Cantera')
        rmpath(dirs{i});
        continue;
    end
    if strfind(dirs{i}, 'cantera')
        rmpath(dirs{i});
    end
end

% Add the Cantera toolbox to the Matlab path
path(path, fullfile(pwd, 'toolbox'));
path(path, fullfile(pwd, 'toolbox', '1D'));

cantera_root = fullfile(pwd, '..', '..');

% Copy the Cantera shared library from the build directory if necessary
if strcmp(getenv('OS'), 'Windows_NT')
    copyfile(fullfile(cantera_root, 'build', 'lib', 'cantera_shared.dll'), ...
        fullfile(pwd, 'toolbox'))
end

% Set path to Python module
if strcmp(getenv('PYTHONPATH'), '')
    setenv('PYTHONPATH', fullfile(cantera_root, 'build', 'python'))
end

% A simple test to make sure that the ctmethods.mex file is present and working
f = Func('polynomial', 3, [1,2,3,4]);
if f(1) == 10
    disp('Cantera MEX file successfully loaded.')
else
    disp('Something is wrong with the Cantera MEX file.')
end

adddir(fullfile(cantera_root, 'data', 'inputs'))
