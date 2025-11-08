classdef AxisymmetricFlow < ct.Flow1D
    % Create an axisymmetric flow domain. ::
    %
    %     >> m = ct.AxisymmetricFlow(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String, ID of the flow.

    methods

        function m = AxisymmetricFlow(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "axisymmetric-flow"
            end

            m@ct.Flow1D('axisymmetric-flow', phase, name);

        end

    end

end
