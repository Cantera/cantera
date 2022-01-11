function v = canteraGitCommit()
    % Get Cantera Git commit hash
    checklib;
    buflen = calllib(ct, 'ct_getGitCommit', 0, '');
    aa = char(zeros(1, buflen));
    [~, aa] = calllib(ct, 'ct_getGitCommit', buflen, aa);
    v = aa;
end
