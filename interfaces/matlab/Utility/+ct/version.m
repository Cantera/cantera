function v = version()
    % Get Cantera version information. ::
    %
    %     >> ct.version()
    %
    % :return:
    %     A string containing the Cantera version.

    v = ct.impl.getString('mCt_version');
end
