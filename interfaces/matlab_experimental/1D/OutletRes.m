function m = OutletRes(id)
    % Create an outlet reservoir domain.
    % m = OutletRes(id)
    %
    % :return:
    %     Instance of :mat:func:`Domain1D` representing an outlet
    %     reservoir.
    %
    m = Domain1D('OutletRes');

    if nargin == 0
        m.setID('outletres');
    else
        m.setID(id);
    end

end
