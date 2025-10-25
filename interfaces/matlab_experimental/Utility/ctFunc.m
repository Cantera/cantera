function output = ctFunc(funcName, varargin)
    % Calls Cantera library functions with single outputs and returns
    % errors if necessary.

    output = clib.ctMatlab.(funcName)(varargin{:});
    % Convert output to double for better data compatibility
    output = double(output);

    if ismember(output, ctErrorCode)
        error('Cantera:ctError', ctGetErr);
    end

end
