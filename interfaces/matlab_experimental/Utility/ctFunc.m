function output = ctFunc(varargin)
    % Calls Cantera library functions with single outputs and returns
    % errors if necessary.

    errorcode = [-1, -999.999, double(intmax('uint64'))];

    funcName = varargin{1};
    output = calllib(ctLib, funcName, varargin{2:end});

    if ismember(output, errorcode)
        error('Cantera:ctError', ctGetErr);
    end

end
