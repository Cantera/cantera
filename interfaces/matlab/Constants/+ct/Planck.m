function r = Planck
    % Get Planck's constant in J·s. ::
    %
    %     >> r = ct.Planck
    %
    % :return:
    %     Planck's constant in J·s.

    r = ct.impl.call('mCt_Planck');
end
