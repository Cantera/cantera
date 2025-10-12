classdef UnstrainedFlow < Flow1D
    % Create an unstrained flow domain. ::
    %
    %     >> m = UnstrainedFlow(phase, id)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`UnstrainedFlow`.

    methods

        function m = UnstrainedFlow(phase, id)
            % Constructor

            if nargin < 2
                id = 'unstrained-flow';
            end

            m@Flow1D('unstrained-flow', phase, id);

        end

    end

end
