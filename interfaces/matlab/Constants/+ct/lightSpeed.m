function r = lightSpeed
    % Get the speed of light in m/s. ::
    %
    %     >> r = ct.lightSpeed
    %
    % :return:
    %     The speed of light in m/s.

    r = ct.impl.call('mCt_lightSpeed');
end
