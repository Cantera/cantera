function p = OneAtm
    % Get one atmosphere in Pa. ::
    %
    %     >> p = ct.OneAtm
    %
    % :return:
    %     One atmosphere in Pascals.

    p = ct.impl.call('mCt_OneAtm');
end
