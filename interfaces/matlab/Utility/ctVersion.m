function v = ctVersion()
    % Get Cantera version information. ::
    %
    %     >> ctVersion()
    %
    % :return:
    %     A string containing the Cantera version.

    v = ctString('mCt_version');
end
