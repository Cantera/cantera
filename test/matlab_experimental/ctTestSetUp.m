function ctTestSetUp()
    % ctTestSetUp
    % Configure paths for running the Cantera MATLAB unit tests and
    % load the Cantera C library.

    thisFile = mfilename('fullpath');
    canteraRoot = fileparts(fileparts(fileparts(thisFile)));

    addpath(genpath(fullfile(canteraRoot, 'interfaces', 'matlab_experimental')));
    addpath(genpath(fullfile(canteraRoot, 'test', 'data')));
    addpath(genpath(fullfile(canteraRoot, 'test', 'matlab_experimental')));

    ctLoad();
end
