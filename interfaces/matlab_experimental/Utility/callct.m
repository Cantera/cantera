function output = callct(varargin)
    % Calls Cantera library functions with single outputs and returns
    % errors if necessary.

    errorcode = [-1, -999.999, double(intmax('uint64'))];

    funcName = varargin{1};
    output = calllib(ct, funcName, varargin{2:end});

    if ismember(output, errorcode)
        error(geterr);
    end

end
