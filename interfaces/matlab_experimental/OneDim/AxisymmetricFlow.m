classdef AxisymmetricFlow < Flow1D
    % Create an axisymmetric flow domain. ::
    %
    %     >> m = AxisymmetricFlow(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param name:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`AxisymmetricFlow`.

    methods

        function m = AxisymmetricFlow(phase, name)
            % Constructor
            arguments
                phase (1,1) Solution
                name (1,1) string = "axisymmetric-flow"
            end

            m@Flow1D('axisymmetric-flow', phase, name);

        end

    end

end
