function m = Outlet(id)
    % Create an outlet domain.
    % m = Outlet(id)
    %
    % :param id:
    %     String ID of the outlet.
    % :return:
    %     Instance of :mat:func:`Domain1D` representing an outlet.
    %
    m = Domain1D('Outlet1D');
    if nargin == 0
        m.setID('outlet');
    else
        m.setID(id);
    end
end
