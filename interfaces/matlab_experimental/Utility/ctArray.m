function output = ctArray(funcName, varargin)
    % Calls Cantera library functions with single outputs and returns
    % errors if necessary.

    errorcode = [-1, -999.999, double(intmax('uint64'))];

    [iok, arr] = clib.ctMatlab.(funcName)(varargin{:});

    if ismember(iok, errorcode)
        error('Cantera:ctError', ctGetErr);
    end

    % Convert output to double for better data compatibility
    output = double(arr);
end
