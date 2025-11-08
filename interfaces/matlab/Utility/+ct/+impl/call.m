function output = call(funcName, varargin)
    % Calls Cantera library functions with single outputs and returns
    % errors if necessary.

    output = clib.ctMatlab.(funcName)(varargin{:});
    % Convert output to double for better data compatibility
    output = double(output);

    if ismember(output, ct.impl.errorCode)
        error('Cantera:ctError', ct.impl.getError());
    end

end
