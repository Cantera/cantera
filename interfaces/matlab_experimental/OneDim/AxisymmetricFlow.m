classdef AxisymmetricFlow < Flow1D
    % Create an axisymmetric flow domain. ::
    %
    %     >> m = AxisymmetricFlow(phase, id)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`AxisymmetricFlow`.

    methods

        function m = AxisymmetricFlow(phase, id)
            % Constructor

            if nargin < 2
                id = 'axisymmetric-flow';
            end

            m@Flow1D('axisymmetric-flow', phase, id);

        end

    end

end
