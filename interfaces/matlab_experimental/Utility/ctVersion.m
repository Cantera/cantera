function v = ctVersion()
    % Get Cantera version information. ::
    %
    %     >> ctVersion()
    %
    % :return:
    %     A string containing the Cantera version.

    ctIsLoaded;
    v = ctString('ct_getCanteraVersion');
end
