function d = getDataDirectories()
    % Get a cell array of the directories searched for data files.
    checklib;
    buflen = calllib(ct, 'ct_getDataDirectories', 0, '', ';');
    aa = char(zeros(1, buflen));
    [~, aa, ~] = calllib(ct, 'ct_getDataDirectories', buflen, aa, ';');
    d = aa;
end
