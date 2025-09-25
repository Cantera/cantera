function ctTestTearDown()
    % ctTestTearDown
    % Unload the Cantera C library and remove test paths from MATLAB.

    ctUnload();

    thisFile = mfilename('fullpath');
    canteraRoot = fileparts(fileparts(fileparts(thisFile)));

    rmpath(genpath(fullfile(canteraRoot, 'interfaces', 'matlab_experimental')));
    rmpath(genpath(fullfile(canteraRoot, 'test', 'data')));
    rmpath(genpath(fullfile(canteraRoot, 'test', 'matlab_experimental')));
end
