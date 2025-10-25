function r = FaradayConstant
    % Get the Faraday constant in C/kmol of electron. ::
    %
    %     >> r = FaradayConstant
    %
    % :return:
    %     The Faraday constant in C/kmol of electron.

    r = ctFunc('ct_Faraday');
end
