classdef Outlet1D < ct.Boundary1D
    % Create an outlet domain. ::
    %
    %     >> m = Outlet1D(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String ID of the outlet.

    methods

        function m = Outlet1D(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "outlet"
            end

            m@ct.Boundary1D('outlet', phase, name);
        end

    end

end
