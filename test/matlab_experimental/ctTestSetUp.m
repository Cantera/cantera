function ctTestSetUp()
    % ctTestSetUp
    % Configure paths for running the Cantera MATLAB unit tests and
    % load the Cantera C library.

    thisFile = mfilename('fullpath');
    canteraRoot = fileparts(fileparts(fileparts(thisFile)));

    testPaths.libPath     = fullfile(canteraRoot, 'build', 'lib');
    testPaths.includePath = fullfile(canteraRoot, 'include');
    testPaths.toolboxPath = fullfile(canteraRoot, 'interfaces', 'matlab_experimental');

    addpath(genpath(testPaths.toolboxPath));
    addpath(genpath(fullfile(canteraRoot, 'test')));

    ctLoad(testPaths);
end
