function r = permeability0
    % Get the vacuum permeability in N/A². ::
    %
    %     >> r = ct.permeability0
    %
    % :return:
    %     The vacuum permeability in N/A².

    r = ct.impl.call('mCt_permeability0');
end
