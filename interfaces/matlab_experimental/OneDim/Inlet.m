classdef Inlet < Boundary1D
    % Create an inlet domain. ::
    %
    %     >> m = Inlet(phase, id)
    %
    % Note that an inlet can only be a terminal domain - it must be
    % either the leftmost or rightmost domain in a stack.
    %
    % :param phase:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String name of the inlet.
    % :return:
    %     Instance of class :mat:class:`Inlet`.

    methods

        function m = Inlet(phase, id)
            % Constructor
            if nargin < 2
                id = 'inlet';
            end

            m@Boundary1D('inlet', phase, id);

        end

    end

end
