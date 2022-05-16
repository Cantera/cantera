function m = Outlet(id)
    % Create an outlet domain.
    m = Domain1D('Outlet1D');
    if nargin == 0
        m.setID('outlet');
    else
        m.setID(id);
    end
end
