function output = ctFunc(funcName, varargin)
    % Calls Cantera library functions with single outputs and returns
    % errors if necessary.

    errorcode = [-1, -999, -999.999, double(intmax('uint64'))];

    output = clib.ctMatlab.(funcName)(varargin{:});
    % Convert output to double for better data compatibility
    output = double(output);

    if ismember(output, errorcode)
        error('Cantera:ctError', ctGetErr);
    end

end
