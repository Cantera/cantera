classdef Outlet < Boundary1D
    % Create an outlet domain. ::
    %
    %     >> m = Outlet(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String ID of the outlet.

    methods

        function m = Outlet(phase, name)
            arguments
                phase (1,1) Solution
                name (1,1) string = "outlet"
            end

            m@Boundary1D('outlet', phase, name);
        end

    end

end
