classdef OutletRes < Boundary1D
    % Create an outlet reservoir domain. ::
    %
    %     >> m = OutletRes(phase, id)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String ID of the outlet reservoir.
    % :return:
    %     Instance of :mat:class:`OutletRes`.

    methods

        function m = OutletRes(phase, id)
            % Constructor

            if nargin < 2
                id = 'outlet-reservoir';
            end

            m@Boundary1D('outlet-reservoir', phase, id);

        end

    end

end
