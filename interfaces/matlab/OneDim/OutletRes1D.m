classdef OutletRes1D < Boundary1D
    % Create an outlet reservoir domain. ::
    %
    %     >> m = OutletRes1D(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String ID of the outlet reservoir.

    methods

        function m = OutletRes1D(phase, name)
            arguments
                phase (1,1) Solution
                name (1,1) string = "outlet-reservoir"
            end

            m@Boundary1D('outlet-reservoir', phase, name);

        end

    end

end
