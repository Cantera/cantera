function v = canteraVersion()
    % Get Cantera version information. ::
    %
    %     >> canteraVersion()
    %
    % :return:
    %     A string containing the Cantera version.
    %
    checklib;
    v = callct2('ct_getCanteraVersion');
end
