classdef FreeFlow < ct.Flow1D
    % Create a free flow domain. ::
    %
    %     >> m = ct.FreeFlow(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String, ID of the flow.

    methods

        function obj = FreeFlow(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "free-flow"
            end

            obj@ct.Flow1D('free-flow', phase, name);

        end

    end

end
