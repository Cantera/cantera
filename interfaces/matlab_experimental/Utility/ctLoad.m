function ctLoad(mode)
    % Load the MATLAB-to-C++ Interface
    %
    % Usage:
    %   ctLoad                % defaults to 'outofprocess'
    %   ctLoad('inprocess')   % load in-process
    %   ctLoad('outofprocess')% load out-of-process

    arguments
        mode (1,1) string {mustBeMember(mode,["inprocess","outofprocess"])} = "outofprocess"
    end

    global ct
    if isempty(ct)
        ct = clibConfiguration("ctMatlab", ExecutionMode=mode);
    end

    fprintf('Cantera %s is ready for use (%s mode).\n', ctVersion, mode);
end
