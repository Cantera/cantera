classdef Outlet < Boundary1D
    % Create an outlet domain. ::
    %
    %     >> m = Outlet(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String ID of the outlet.
    % :return:
    %     Instance of :mat:class:`Outlet`.

    methods

        function m = Outlet(phase, name)
            % Constructor
            arguments
                phase (1,1) Solution
                name (1,1) string = "outlet"
            end

            m@Boundary1D('outlet', phase, name);
        end

    end

end
