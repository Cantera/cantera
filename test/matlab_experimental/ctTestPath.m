% ctTestPath.m
% Set up environment for testing the Cantera Matlab interface
% from within the Cantera source tree. Run this file from the
% root of the Cantera source tree, for example:
%
%    cd ~/src/cantera
%    run interfaces/matlab_experimental/testpath.m

% get the list of directories on the Matlab path
dirs = split(path, pathsep);

% if 'cantera' is already in the path, remove it
for i = 1:length(dirs)
    if contains(dirs{i}, 'CANTERA', 'IgnoreCase', true)
        rmpath(dirs{i});
    end
end

cantera_root = getenv('CANTERA_ROOT');

% Add the Cantera toolbox to the Matlab path
addpath(genpath([cantera_root, '/interfaces/matlab_experimental']));
addpath(genpath([cantera_root, '/test/matlab_experimental']));
addpath(genpath([cantera_root, '/test/data']));
