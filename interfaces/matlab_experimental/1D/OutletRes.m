function m = OutletRes(id)
    % Create an outlet reservoir domain.
    m = Domain1D(-2);
    if nargin == 0
        m.setID('outletres');
    else
        m.setID(id);
    end
end
