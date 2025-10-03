function v = ctVersion()
    % Get Cantera version information. ::
    %
    %     >> ctVersion()
    %
    % :return:
    %     A string containing the Cantera version.

    v = ctString('ct_version');
end
