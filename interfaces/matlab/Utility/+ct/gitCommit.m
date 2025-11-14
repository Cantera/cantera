function v = gitCommit()
    % Get Cantera Git commit hash. ::
    %
    %     >> ct.gitCommit()
    %
    % :return:
    %     A string containing the Git commit hash for the current version of Cantera.

    ct.isLoaded(true);
    v = ct.impl.getString('mCt_gitCommit');
end
