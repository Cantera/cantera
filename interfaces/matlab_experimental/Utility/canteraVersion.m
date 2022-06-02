function v = canteraVersion()
    % Get Cantera version information
    % canteraVersion()
    %
    % :return:
    %     A string containing the Cantera version
    %
    checklib;
    buflen = calllib(ct, 'ct_getCanteraVersion', 0, '');
    aa = char(zeros(1, buflen));
    [~, aa] = calllib(ct, 'ct_getCanteraVersion', buflen, aa);
    v = aa;
end
