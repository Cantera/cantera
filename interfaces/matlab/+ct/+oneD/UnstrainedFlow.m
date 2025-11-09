classdef UnstrainedFlow < ct.oneD.Flow
    % Create an unstrained flow domain. ::
    %
    %     >> m = ct.oneD.UnstrainedFlow(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`ct.oneD.UnstrainedFlow`.

    methods

        function obj = UnstrainedFlow(phase, name)
            % Constructor
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "unstrained-flow"
            end

            obj@ct.oneD.Flow('unstrained-flow', phase, name);

        end

    end

end
