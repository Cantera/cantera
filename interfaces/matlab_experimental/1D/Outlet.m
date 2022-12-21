classdef Outlet < Domain1D
    % Create an outlet domain.
    %
    % m = Outlet(id)
    %
    % :param id:
    %     String ID of the outlet.
    % :return:
    %     Instance of :mat:class:`Outlet`.
    %

    methods

        % Constructor
        function m = Outlet(id)

            m = m@Domain1D('Outlet1D');

            if nargin == 0
                m.setID('outlet');
            else
                m.setID(id);
            end

        end

    end

end
