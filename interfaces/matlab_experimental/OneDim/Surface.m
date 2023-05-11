classdef Surface < Boundary1D
    % Create a surface domain. ::
    %
    %     >> m = Surface(phase, id)
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String ID of surface
    % :return:
    %     Instance of class :mat:class:`Surface` representing a
    %     non-reacting surface.

    methods

        %% Surface Class Constructor

        function m = Surface(phase, id)

            if nargin < 2
                id= 'surface';
            end

            m@Boundary1D('surface', phase, id);

        end

    end

end
