function r = ElectronMass
    % Get the electron mass in kg. ::
    %
    %     >> r = ct.ElectronMass
    %
    % :return:
    %     The electron mass in kg.

    r = ct.impl.call('mCt_ElectronMass');
end
