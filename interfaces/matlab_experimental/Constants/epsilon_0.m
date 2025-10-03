function r = epsilon_0
    % Get the vacuum permittivity in C^2/N/m^2. ::
    %
    %     >> r = epsilon_0
    %
    % :return:
    %     The vacuum permittivity in C^2/N/m^2.

    r = ctFunc('ct_epsilon_0');
end
