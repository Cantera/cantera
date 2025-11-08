function r = GasConstant
    % Get the universal gas constant in J/kmol/K. ::
    %
    %     >> r = gasConstant
    %
    % :return:
    %     The universal gas constant in J/kmol/K.

    r = ctFunc('mCt_GasConstant');
end
