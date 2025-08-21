function ctTestTearDown()
    % ctTestTearDown
    % Unload the Cantera C library and remove test paths from MATLAB.

    ctUnload();

    thisFile = mfilename('fullpath');
    canteraRoot = fileparts(fileparts(fileparts(thisFile)));

    toolboxPath = fullfile(canteraRoot, 'interfaces', 'matlab_experimental');
    testPath    = fullfile(canteraRoot, 'test');

    rmpath(genpath(toolboxPath));
    rmpath(genpath(testPath));
end
