function v = ctGitCommit()
    % Get Cantera Git commit hash. ::
    %
    %     >> ctGitCommit()
    %
    % :return:
    %     A string containing the Git commit hash for the current version of Cantera.

    ctIsLoaded;
    v = ctString('ct_getGitCommit');
end
