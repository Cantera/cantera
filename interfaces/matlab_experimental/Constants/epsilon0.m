function r = epsilon0
    % Get the vacuum permittivity in C^2/N/m^2. ::
    %
    %     >> r = epsilon0
    %
    % :return:
    %     The vacuum permittivity in C^2/N/m^2.

    r = ctFunc('ct_epsilon0');
end
