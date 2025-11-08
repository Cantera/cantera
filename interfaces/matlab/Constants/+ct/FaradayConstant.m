function r = FaradayConstant
    % Get the Faraday constant in C/kmol of electron. ::
    %
    %     >> r = ct.FaradayConstant
    %
    % :return:
    %     The Faraday constant in C/kmol of electron.

    r = ct.impl.call('mCt_Faraday');
end
