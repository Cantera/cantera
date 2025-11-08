function r = Boltzmann
    % Get Boltzmann's constant in J/K. ::
    %
    %     >> r = Boltzmann
    %
    % :return:
    %     Boltzmann's constant in J/K.

    r = ct.impl.call('mCt_Boltzmann');
end
