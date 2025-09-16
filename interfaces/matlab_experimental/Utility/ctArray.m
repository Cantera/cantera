function output = ctArray(funcName, buflen, varargin)
    % Calls Cantera library functions with single outputs and returns
    % errors if necessary.

    errorcode = [-1, -999.999, double(intmax('uint64'))];

    buf = clib.array.ctMatlab.Double(buflen);
    iok = clib.ctMatlab.(funcName)(varargin{:}, buf);

    iok = double(iok);
    if ismember(iok, errorcode)
        error('Cantera:ctError', ctGetErr);
    end

    % Convert output to double for better data compatibility
    output = buf.double;
end
