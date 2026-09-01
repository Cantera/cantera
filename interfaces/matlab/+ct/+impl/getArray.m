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

    persistent arrayCache
    if isempty(arrayCache)
        arrayCache = containers.Map('KeyType','double','ValueType','any');
    end

    % Get or create buffer for this size
    if arrayCache.isKey(buflen)
        buf = arrayCache(buflen);
    else
        buf = clib.array.ctMatlab.Double(buflen);
        arrayCache(buflen) = buf;
    end

    iok = clib.ctMatlab.(funcName)(args{:}, buf, extraArgs{:});

    ct.impl.checkErrorCode(iok);

    % Convert output to double for better data compatibility
    output = buf.double;
end
