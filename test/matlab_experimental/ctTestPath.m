% ctTestPath.m
% Set up environment for testing the Cantera Matlab interface
% from within the Cantera source tree. Run this file from the
% root of the Cantera source tree, for example:
%
%    cd ~/src/cantera
%    run interfaces/matlab_experimental/testpath.m

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

cantera_root = fullfile(pwd);

% Copy the Cantera shared library from the build directory if necessary
if ispc
    libname = 'cantera_shared.dll';
elseif ismac
    libname = 'libcantera_shared.dylib';
elseif isunix
    libname = 'libcantera_shared.so';
end

copyfile(fullfile(cantera_root, 'build', 'lib', libname), ...
         fullfile(pwd, 'test', 'matlab_experimental'));

% Add the Cantera toolbox to the Matlab path
addpath(genpath([cantera_root, '/interfaces/matlab_experimental']));
addpath(genpath([cantera_root, '/test/matlab_experimental']));

% Set path to Python module
if strcmp(getenv('PYTHONPATH'), '')
    setenv('PYTHONPATH', fullfile(cantera_root, 'build', 'python'))
end
