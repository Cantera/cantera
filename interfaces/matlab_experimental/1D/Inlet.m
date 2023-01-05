classdef Inlet < Domain1D
    % Create an inlet domain. ::
    %
    %     >> m = Inlet(id)
    %
    % Note that an inlet can only be a terminal domain - it must be
    % either the leftmost or rightmost domain in a stack.
    %
    % :param id:
    %     String name of the inlet.
    % :return:
    %     Instance of class :mat:class:`Inlet`.
    %

    methods

        function m = Inlet(id)
            % Constructor

            m@Domain1D('Inlet1D');

            if nargin == 0
                m.setID('inlet');
            else
                m.setID(id);
            end

        end

    end

end
