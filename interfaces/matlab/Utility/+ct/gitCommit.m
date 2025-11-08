function v = gitCommit()
    % Get Cantera Git commit hash. ::
    %
    %     >> ct.gitCommit()
    %
    % :return:
    %     A string containing the Git commit hash for the current version of Cantera.

    isLoaded(true);
    v = ctString('mCt_getGitCommit');
end
