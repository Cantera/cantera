function r = StefanBoltz
    % Get Stefan-Boltzmann constant in W/m²/K⁴. ::
    %
    %     >> r = StefanBoltz
    %
    % :return:
    %     Stefan-Boltzmann constant in W/m²/K⁴.

    r = ct.impl.call('mCt_StefanBoltz');
end
