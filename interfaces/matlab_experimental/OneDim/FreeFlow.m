classdef FreeFlow < Flow1D
    % Create a free flow domain. ::
    %
    %     >> m = FreeFlow(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`FreeFlow`.

    methods

        function m = FreeFlow(phase, name)
            % Constructor
            arguments
                phase (1,1) Solution
                name (1,1) string = "free-flow"
            end

            m@Flow1D('free-flow', phase, name);

        end

    end

end
