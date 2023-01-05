classdef OutletRes < Domain1D
    % Create an outlet reservoir domain. ::
    %
    %     >> m = OutletRes(id)
    %
    % :param id:
    %     String ID of the outlet reservoir.
    % :return:
    %     Instance of :mat:class:`OutletRes`.
    %

    methods

        function m = OutletRes(id)
            % Constructor

            m@Domain1D('OutletRes');

            if nargin == 0
                m.setID('outletres');
            else
                m.setID(id);
            end

        end

    end

end
