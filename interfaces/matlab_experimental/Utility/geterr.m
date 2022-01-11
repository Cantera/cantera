function e = geterr()
    % Get the error message from Cantera error.
    checklib;
    try
        buflen = calllib(ct, 'ct_getCanteraError', 0, '');
        aa = char(zeros(1, buflen));
        [~, aa] = calllib(ct, 'ct_getCanteraError', buflen, aa);
        e = aa;
    catch ME
        e = getReport(ME);
    end
end
