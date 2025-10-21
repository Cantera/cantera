function ctLoad(mode)
    % Load the MATLAB-to-C++ Interface
    %
    % Usage:
    %   ctLoad                % defaults to 'outofprocess'
    %   ctLoad('inprocess')   % load in-process
    %   ctLoad('outofprocess')% load out-of-process

    arguments
        mode (1,1) string {mustBeMember(mode, ...
                           ["inprocess","outofprocess"])} = "outofprocess"
    end

    pathVar = dictionary();
    if ismac
        arch = computer("arch");
        matlabLibPath = matlabroot + "/bin/" + arch + pathsep + ...
                        matlabroot + "/sys/os/" + arch;
        pathVar = dictionary("DYLD_LIBRARY_PATH", matlabLibPath);
    end

    global ct
    if isempty(ct)
        if isMATLABReleaseOlderThan("R2025a") || strcmp(mode, "inprocess")
            ct = clibConfiguration("ctMatlab", ExecutionMode=mode);
        else
            ct = clibConfiguration("ctMatlab", ExecutionMode=mode, ...
                                   OutOfProcessEnvironmentVariables=pathVar);
        end
    end

    fprintf('Cantera %s is ready for use (%s mode).\n', ctVersion, mode);
end
