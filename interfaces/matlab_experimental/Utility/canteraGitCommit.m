function v = canteraGitCommit()
    % Get Cantera Git commit hash
    % canteraGitCommit()
    %
    % :return:
    %     A string containing the Git commit hash for the current version of Cantera
    %
    checklib;
    buflen = callct('ct_getGitCommit', 0, '');
    aa = char(zeros(1, buflen));
    [~, aa] = callct('ct_getGitCommit', buflen, aa);
    v = aa;
end
