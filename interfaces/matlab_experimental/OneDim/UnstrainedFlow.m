classdef UnstrainedFlow < Flow1D
    % Create an unstrained flow domain. ::
    %
    %     >> m = UnstrainedFlow(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`UnstrainedFlow`.

    methods

        function m = UnstrainedFlow(phase, name)
            % Constructor
            arguments
                phase (1,1) Solution
                name (1,1) string = "unstrained-flow"
            end

            m@Flow1D('unstrained-flow', phase, name);

        end

    end

end
