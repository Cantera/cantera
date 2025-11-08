function r = Avogadro
    % Get Avogadro's number in 1/kmol. ::
    %
    %     >> r = ct.Avogadro
    %
    % :return:
    %     Avogadro's number in 1/kmol.

    r = ct.impl.call('mCt_Avogadro');
end
