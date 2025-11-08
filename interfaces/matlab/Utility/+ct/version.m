function v = version()
    % Get Cantera version information. ::
    %
    %     >> ct.version()
    %
    % :return:
    %     A string containing the Cantera version.

    v = ctString('mCt_version');
end
