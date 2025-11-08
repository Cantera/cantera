function r = epsilon0
    % Get the vacuum permittivity in C²/N/m². ::
    %
    %     >> r = ct.epsilon0
    %
    % :return:
    %     The vacuum permittivity in C²/N/m².

    r = ct.impl.call('mCt_epsilon0');
end
