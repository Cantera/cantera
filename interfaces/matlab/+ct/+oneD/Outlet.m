classdef Outlet < ct.oneD.Boundary
    % Create an outlet domain. ::
    %
    %     >> m = ct.oneD.Outlet(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String ID of the outlet.

    methods

        function obj = Outlet(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "outlet"
            end

            obj@ct.oneD.Boundary('outlet', phase, name);
        end

    end

end
