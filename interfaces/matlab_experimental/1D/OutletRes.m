function m = OutletRes(id)
    % Create an outlet reservoir domain.
    m = Domain1D('OutletRes');
    if nargin == 0
        m.setID('outletres');
    else
        m.setID(id);
    end
end
