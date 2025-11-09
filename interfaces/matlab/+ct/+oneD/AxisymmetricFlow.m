classdef AxisymmetricFlow < ct.oneD.Flow
    % Create an axisymmetric flow domain. ::
    %
    %     >> m = ct.oneD.AxisymmetricFlow(phase, name)
    %
    % :param phase:
    %     Instance of class :mat:class:`ct.Solution`.
    % :param name:
    %     String, ID of the flow.

    methods

        function obj = AxisymmetricFlow(phase, name)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "axisymmetric-flow"
            end

            obj@ct.oneD.Flow('axisymmetric-flow', phase, name);

        end

    end

end
