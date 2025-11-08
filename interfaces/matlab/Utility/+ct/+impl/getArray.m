function output = getArray(funcName, buflen, varargin)
    % Calls Cantera library functions that return arrays and handles
    % errors if necessary.

    idx = find(strcmp(varargin, 'extraArgs'), 1);

    if ~isempty(idx)
        args = varargin(1:idx-1);
        extraArgs = varargin(idx+1:end);
    else
        args = varargin;
        extraArgs = {};
    end

    buf = clib.array.ctMatlab.Double(buflen);
    iok = clib.ctMatlab.(funcName)(args{:}, buf, extraArgs{:});

    iok = double(iok);
    if ismember(iok, ct.impl.errorCode)
        error('Cantera:ctError', ct.impl.getError());
    end

    % Convert output to double for better data compatibility
    output = buf.double;
end
