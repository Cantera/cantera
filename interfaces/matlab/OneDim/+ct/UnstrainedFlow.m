classdef UnstrainedFlow < ct.Flow1D
    % Create an unstrained flow domain. ::
    %
    %     >> m = ct.UnstrainedFlow(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`ct.UnstrainedFlow`.

    methods

        function m = UnstrainedFlow(phase, name)
            % Constructor
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "unstrained-flow"
            end

            m@ct.Flow1D('unstrained-flow', phase, name);

        end

    end

end
