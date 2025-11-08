function r = Boltzmann
    % Get Boltzmann's constant in J/K. ::
    %
    %     >> r = ct.Boltzmann
    %
    % :return:
    %     Boltzmann's constant in J/K.

    r = ct.impl.call('mCt_Boltzmann');
end
