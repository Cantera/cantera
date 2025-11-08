function load(mode)
    % Load the MATLAB-to-C++ Interface
    %
    % Usage:
    %   ct.load                % defaults to 'outofprocess'
    %   ct.load('inprocess')   % load in-process
    %   ct.load('outofprocess')% load out-of-process

    arguments
        mode (1,1) string {mustBeMember(mode, ...
                           ["inprocess", "outofprocess"])} = "outofprocess"
    end

    if isLoaded
        mode_ = executionMode;
        msg = sprintf('Cantera %s is already loaded (%s mode).', version, mode_);
        warning("load:IsLoaded", msg);
        if mode == mode_ && mode == "inprocess"
            return
        elseif mode ~= mode_ && mode_ == "inprocess"
            error("load:LoadFailed", ...
                  ("Unloading of `ctMatlab` library is not supported for " + ...
                   "'inprocess' execution mode. Restart MATLAB to update."));
            return
        else
            ctUnload;
        end
    end

    pathVar = dictionary();
    if ismac
        arch = computer("arch");
        matlabLibPath = matlabroot + "/bin/" + arch + pathsep + ...
                        matlabroot + "/sys/os/" + arch;
        pathVar = dictionary("DYLD_LIBRARY_PATH", matlabLibPath);
    end

    global ctMatlab
    if ~isLoaded
        if isMATLABReleaseOlderThan("R2025a") || strcmp(mode, "inprocess")
            ctMatlab = clibConfiguration("ctMatlab", ExecutionMode=mode);
        else
            ctMatlab = clibConfiguration("ctMatlab", ExecutionMode=mode, ...
                                   OutOfProcessEnvironmentVariables=pathVar);
        end
    end

    fprintf('Cantera %s is ready for use (%s mode).\n', version, mode);
end
