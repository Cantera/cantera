function p = OneAtm
    % Get one atmosphere in Pa. ::
    %
    %     >> p = oneatm
    %
    % :return:
    %     One atmosphere in Pascals.

    p = ct.impl.call('mCt_OneAtm');
end
