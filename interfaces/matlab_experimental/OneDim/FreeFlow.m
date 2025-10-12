classdef FreeFlow < Flow1D
    % Create a free flow domain. ::
    %
    %     >> m = FreeFlow(phase, id)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`FreeFlow`.

    methods

        function m = FreeFlow(phase, id)
            % Constructor

            if nargin < 2
                id = 'free-flow';
            end

            m@Flow1D('free-flow', phase, id);

        end

    end

end
