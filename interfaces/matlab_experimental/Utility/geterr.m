function e = geterr()
    % Get the error message from a Cantera error.
    %
    checklib;
    try
        buflen = calllib(ct, 'ct_getCanteraError', 0, '');
        aa = char([zeros(1, buflen, 'int8')]);
        ptr = libpointer('cstring', aa);
        [~, bb] = calllib(ct, 'ct_getCanteraError', buflen, ptr);
        e = bb;
        clear aa bb ptr
    catch ME
        e = getReport(ME);
    end
end
