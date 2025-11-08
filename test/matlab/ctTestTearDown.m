function ctTestTearDown()
    % ctTestTearDown
    % Unload the Cantera C library and remove test paths from MATLAB.

    if ct.executionMode == "outofprocess"
        ct.unload;
    else
        ct.cleanUp;
    end

    thisFile = mfilename('fullpath');
    canteraRoot = fileparts(fileparts(fileparts(thisFile)));

    rmpath(genpath(fullfile(canteraRoot, 'interfaces', 'matlab')));
    rmpath(genpath(fullfile(canteraRoot, 'test', 'data')));
    rmpath(genpath(fullfile(canteraRoot, 'test', 'matlab')));
end
