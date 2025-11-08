function r = ElectronCharge
    % Get the elementary charge in Coulombs. ::
    %
    %     >> r = ct.ElectronCharge
    %
    % :return:
    %     The elementary charge in Coulombs.

    r = ct.impl.call('mCt_ElectronCharge');
end
