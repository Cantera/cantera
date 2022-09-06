function v = canteraGitCommit()
    % Get Cantera Git commit hash
    % canteraGitCommit()
    %
    % :return:
    %     A string containing the Git commit hash for the current version of Cantera
    %
    checklib;
    v = callct2('ct_getGitCommit');
end
